// sr_core.h
// Search core shared by the segment-routing heuristics (amns.cpp and
// classical_solution.cpp). Only the per-slot drivers differ between them.
//
//   Router        : all-pairs shortest distances for one slot topology, and
//                   ECMP splits of a segment (u, v): every node of the
//                   forwarding graph FG(u, v) divides its incoming flow evenly
//                   between its outgoing shortest-path arcs, as the checker does.
//   SlotSearch    : per-slot state (waypoints, fractional link flow, budget,
//                   lexicographic objective) and the neighborhoods.
//   FutureSlots   : loads of every later slot under the current slot-0
//                   waypoints, used by the robust slot-0 phase.
//   RobustSearch  : slot-0 neighborhood scored against all slots at once.
//
// Everything lives in an anonymous namespace: each solver's translation unit
// gets its own copy.

#ifndef SR_CORE_H
#define SR_CORE_H

#include "loader.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;

namespace {

// Infinity used for "unreachable" shortest paths.
constexpr double INF = numeric_limits<double>::infinity();
constexpr double kEps = 1e-12;
constexpr int kShortlist = 16;       // residual-ranked nodes delta-scored per position
constexpr int kShortlistSlot0 = 32;  // slot 0 is searched more aggressively
constexpr int kDetourTry = 3;        // try_move the best delta-scored candidates
constexpr int kDetourTrySlot0 = 5;
constexpr int kDetourRank = 12;      // random pick range when shaking
constexpr double kCritRel = 0.002;   // near-tie band defining the critical cut
constexpr double kShakeSlack = 0.02; // shake may raise max_load by up to 2%
constexpr double kSpreadRel = 0.5;   // links above this fraction of max_load are "warm"
constexpr int kExtendDemands = 48;   // demands examined per extend/robust call (rotating window)
constexpr int kRobustWide = 64;      // detour nodes scored per position for demands on the all-slot peak
constexpr int kRobustPair = 16;      // partner nodes scored per anchor in a two-waypoint detour
constexpr int kRobustKicks = 30;     // robust phase stops after this many kicks without a new best
constexpr int kPairFirst = 6;        // insert2: first-waypoint nodes per position (besides anchors)
constexpr int kPairSecond = 6;       // insert2: partner nodes per first waypoint and side
constexpr size_t kSegCacheCap = 6000000;    // cached split entries per Router
constexpr size_t kFutureCacheCap = 6000000; // cached entries in FutureSlots

// One directed edge in the routing graph.
//   to      : head node (tail node in Router::radj)
//   link_id : original network link id (used for capacity, forbidden set)
//   metric  : routing cost (e.g. latency)
struct Arc {
    int to;
    int link_id;
    double metric;
};

// ECMP split of one segment: (link id, fraction of the segment's flow).
// Each link appears at most once.
using Split = vector<pair<int, double>>;

// Pack an ordered node pair (u, v) into a 64-bit key for hashing.
static uint64_t pack_seg(int u, int v)
{
    return (static_cast<uint64_t>(static_cast<uint32_t>(u)) << 32)
           | static_cast<uint32_t>(v);
}

// True if an arc of cost `metric` into a node at distance `d_head` (to the
// segment end, or from the segment start) lies on a shortest path to/from a
// node at distance `d_tail`. Metrics are integers, so sums are exact; the
// tolerance only guards against rounding on huge totals.
static bool on_sp(double d_head, double metric, double d_tail)
{
    return d_head < INF && fabs(d_head + metric - d_tail) <= 1e-12 * max(1.0, d_tail);
}

// Merge repeated links of a concatenated route (overlapping segments).
static void merge_split(Split& s)
{
    if (s.size() < 2) return;
    sort(s.begin(), s.end());
    size_t k = 0;
    for (size_t i = 1; i < s.size(); ++i) {
        if (s[i].first == s[k].first) s[k].second += s[i].second;
        else s[++k] = s[i];
    }
    s.resize(k + 1);
}

// Scratch buffers for ecmp_split (stamp-marked, never cleared).
struct EcmpScratch {
    vector<double> f;     // flow reaching each node
    vector<double> dv;    // distance to the segment end
    vector<int> mark;
    vector<int> nodes;
    int stamp = 0;

    int next(int n)
    {
        if (static_cast<int>(mark.size()) < n) {
            f.assign(n, 0.0);
            dv.assign(n, INF);
            mark.assign(n, 0);
        }
        if (++stamp == numeric_limits<int>::max()) {
            fill(mark.begin(), mark.end(), 0);
            stamp = 1;
        }
        return stamp;
    }
};

// ECMP split of one unit of flow from u to v. D(x) is the distance from x to
// v in the topology; links with down[link] set are unusable. Nodes of FG(u, v)
// are processed by decreasing distance, so each node's inflow is complete
// before it is divided between its next hops.
template <class DistTo>
static void ecmp_split(int u, int v, const vector<vector<Arc>>& adj,
                       const vector<char>& down, DistTo D,
                       EcmpScratch& sc, Split& out)
{
    out.clear();
    if (u == v) return;
    const int stamp = sc.next(static_cast<int>(adj.size()));
    sc.nodes.clear();
    sc.nodes.push_back(u);
    sc.mark[u] = stamp;
    sc.dv[u] = D(u);
    for (size_t i = 0; i < sc.nodes.size(); ++i) {
        const int x = sc.nodes[i];
        sc.f[x] = 0.0;
        if (x == v) continue;
        for (const Arc& a : adj[x]) {
            if (down[a.link_id] || sc.mark[a.to] == stamp) continue;
            const double dy = D(a.to);
            if (!on_sp(dy, a.metric, sc.dv[x])) continue;
            sc.mark[a.to] = stamp;
            sc.dv[a.to] = dy;
            sc.nodes.push_back(a.to);
        }
    }
    sort(sc.nodes.begin(), sc.nodes.end(),
         [&](int a, int b) { return sc.dv[a] > sc.dv[b]; });
    sc.f[u] = 1.0;
    for (int x : sc.nodes) {
        if (x == v || sc.f[x] <= 0.0) continue;
        int k = 0;
        for (const Arc& a : adj[x])
            if (!down[a.link_id] && sc.mark[a.to] == stamp
                && on_sp(sc.dv[a.to], a.metric, sc.dv[x]))
                ++k;
        if (k == 0) continue;
        const double share = sc.f[x] / k;
        for (const Arc& a : adj[x]) {
            if (down[a.link_id] || sc.mark[a.to] != stamp
                || !on_sp(sc.dv[a.to], a.metric, sc.dv[x]))
                continue;
            out.push_back({a.link_id, share});
            sc.f[a.to] += share;
        }
    }
}

// Precomputes and stores all-pairs shortest distances for a fixed
// "forbidden" set of links (failures for the current time slot), and serves
// cached ECMP splits of segments.
struct Router {
    int num_nodes = 0;         // number of nodes (max_id + 1)
    int max_link_id = -1;      // largest link id seen
    vector<int> node_ids;      // list of valid node ids
    vector<vector<Arc>> adj;   // out-arcs per node
    vector<vector<Arc>> radj;  // in-arcs per node (Arc::to is the tail)
    vector<int> link_from, link_to;  // per link id: endpoints
    vector<char> down;         // per link id: unusable in this slot
    vector<double> dist;       // shortest distance, indexed src*N + v
    vector<int> order;         // per source: reachable nodes by increasing distance
    vector<int> order_len;     // per source: number of reachable nodes

    mutable vector<vector<int>> to_order;          // per target: nodes by increasing distance to it (lazy)
    mutable unordered_map<uint64_t, Split> cache;  // (a, b) -> ECMP split
    mutable size_t cache_size = 0;
    mutable EcmpScratch scratch;

    // Flattened 2D index for (src, v).
    int idx(int src, int v) const { return src * num_nodes + v; }

    // Build the graph from the network topology.
    void build(const Network& net)
    {
        // Nodes are assumed to be sparse ids; allocate arrays up to max id.
        int max_id = -1;
        for (const auto& n : net.nodes)
            if (n.id > max_id) max_id = n.id;
        num_nodes = max_id + 1;

        node_ids.clear();
        node_ids.reserve(net.nodes.size());
        for (const auto& n : net.nodes)
            node_ids.push_back(n.id);

        adj.assign(num_nodes, {});
        radj.assign(num_nodes, {});
        max_link_id = -1;
        for (const auto& l : net.links)
            max_link_id = max(max_link_id, l.id);
        link_from.assign(max_link_id + 1, -1);
        link_to.assign(max_link_id + 1, -1);
        for (const auto& l : net.links) {
            link_from[l.id] = l.from;
            link_to[l.id] = l.to;
            adj[l.from].push_back({l.to, l.id, l.metric});
            radj[l.to].push_back({l.from, l.id, l.metric});
            // Undirected networks get the reverse arc too (same link id).
            if (!net.directed) {
                adj[l.to].push_back({l.from, l.id, l.metric});
                radj[l.from].push_back({l.to, l.id, l.metric});
            }
        }
    }

    // Run Dijkstra from every valid node, skipping forbidden links.
    void compute(const unordered_set<int>& forbidden)
    {
        const int N = num_nodes;
        down.assign(max_link_id + 1, 0);
        for (int lid : forbidden)
            if (lid >= 0 && lid <= max_link_id) down[lid] = 1;
        cache.clear();
        cache_size = 0;
        to_order.assign(N, {});
        dist.assign(static_cast<size_t>(N) * N, INF);
        order.assign(static_cast<size_t>(N) * N, -1);
        order_len.assign(N, 0);

        vector<char> visited(N, 0);
        using P = pair<double, int>; // (distance, node)
        for (int s : node_ids) {
            fill(visited.begin(), visited.end(), 0);
            priority_queue<P, vector<P>, greater<P>> pq;
            dist[idx(s, s)] = 0.0;
            pq.push({0.0, s});

            while (!pq.empty()) {
                auto [d, u] = pq.top();
                pq.pop();
                if (visited[u]) continue; // stale entry
                visited[u] = 1;
                order[idx(s, order_len[s]++)] = u;
                for (const auto& arc : adj[u]) {
                    if (down[arc.link_id]) continue;
                    const double nd = d + arc.metric;
                    const int j = idx(s, arc.to);
                    if (nd < dist[j]) {
                        dist[j] = nd;
                        pq.push({nd, arc.to});
                    }
                }
            }
        }
    }

    // Nodes that reach b, by increasing distance to b.
    const vector<int>& nodes_to(int b) const
    {
        vector<int>& o = to_order[b];
        if (o.empty()) {
            for (int x : node_ids)
                if (dist[idx(x, b)] < INF) o.push_back(x);
            sort(o.begin(), o.end(), [&](int x, int y) {
                return dist[idx(x, b)] < dist[idx(y, b)];
            });
        }
        return o;
    }

    // ECMP split of segment a -> b; nullptr if b is unreachable from a.
    // The pointer stays valid until the next trim_cache() or compute().
    const Split* segment(int a, int b) const
    {
        static const Split kNoLinks;
        if (a == b) return &kNoLinks;
        if (a < 0 || b < 0 || a >= num_nodes || b >= num_nodes)
            return nullptr;
        if (dist[idx(a, b)] == INF)
            return nullptr;
        const uint64_t key = pack_seg(a, b);
        auto it = cache.find(key);
        if (it != cache.end()) return &it->second;
        Split s;
        ecmp_split(a, b, adj, down,
                   [&](int x) { return dist[idx(x, b)]; }, scratch, s);
        cache_size += s.size() + 1;
        return &cache.emplace(key, std::move(s)).first->second;
    }

    // Drop the split cache once it grows too large. Call only between moves,
    // never while holding a pointer from segment().
    void trim_cache() const
    {
        if (cache_size > kSegCacheCap) {
            cache.clear();
            cache_size = 0;
        }
    }

    // Convert a waypoint list [w1, w2, ...] into the ECMP link fractions of
    // the concatenated segments s->w1, w1->w2, ..., wk->t.
    // Returns false if any segment is unreachable.
    bool expand(int s, int t, const vector<int>& waypoints, Split& links) const
    {
        links.clear();
        if (s == t && waypoints.empty())
            return true; // trivial path

        // Build the full node sequence: source, waypoints (dedup consecutive),
        // destination.
        vector<int> seq;
        seq.reserve(waypoints.size() + 2);
        seq.push_back(s);
        for (int w : waypoints) {
            if (w != seq.back())
                seq.push_back(w);
        }
        if (seq.back() != t)
            seq.push_back(t);

        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            const int a = seq[i];
            const int b = seq[i + 1];
            if (a == b) continue;
            const Split* seg = segment(a, b);
            if (!seg) {
                links.clear();
                return false;
            }
            links.insert(links.end(), seg->begin(), seg->end());
        }
        merge_split(links);
        return true;
    }
};

// Sorted, de-duplicated segment pairs (consecutive hops) of a demand route
//   s -> w1 -> w2 -> ... -> t
// Used to measure how many segments changed vs. another route.
static void segment_pairs(int s, int t, const vector<int>& w, vector<uint64_t>& segs)
{
    segs.clear();
    if (s == t && w.empty())
        return;
    int prev = s;
    for (int x : w) {
        if (x == prev) continue; // skip duplicate
        segs.push_back(pack_seg(prev, x));
        prev = x;
    }
    if (prev != t)
        segs.push_back(pack_seg(prev, t));
    sort(segs.begin(), segs.end());
    segs.erase(unique(segs.begin(), segs.end()), segs.end());
}

// Symmetric-difference size between the segment sets of two waypoint lists.
// This is the "number of segments changed" and feeds the budget check.
static int segment_distance(int s, int t,
                            const vector<int>& a,
                            const vector<int>& b)
{
    static vector<uint64_t> sa, sb;   // reused buffers (single-threaded)
    segment_pairs(s, t, a, sa);
    segment_pairs(s, t, b, sb);
    size_t i = 0, j = 0;
    int common = 0;
    while (i < sa.size() && j < sb.size()) {
        if (sa[i] < sb[j]) ++i;
        else if (sb[j] < sa[i]) ++j;
        else { ++common; ++i; ++j; }
    }
    return static_cast<int>(sa.size() + sb.size()) - 2 * common;
}

// Look up the budget value for a given time slot; +inf if slot <= 0
// (i.e. no reconfiguration budget is enforced for the first slot).
static double budget_for_slot(const Instance& inst, int slot)
{
    if (slot <= 0)
        return numeric_limits<double>::infinity();
    for (const auto& b : inst.scenario.budget)
        if (b.t == slot)
            return b.value;
    return 0.0;
}

// Traffic that must cross a link at `slot` whatever the waypoints: from a
// source whose only usable out-link is that link (or into a target whose only
// usable in-link is), followed along single-link chains. A link carrying
// nothing but forced traffic is ignored by the objective (see eff_load):
// no move can change it, and pinning the peak on it would stall the search.
static vector<double> forced_flow_at(const Instance& inst, int num_nodes,
                                     int n_links, int slot)
{
    vector<double> forced(n_links, 0.0);
    unordered_set<int> down;
    for (const auto& iv : inst.scenario.interventions)
        if (iv.t == slot) down.insert(iv.links.begin(), iv.links.end());

    // (link id, other end) of every usable link, per node.
    const int N = num_nodes;
    vector<vector<pair<int, int>>> out(N), in(N);
    for (const auto& l : inst.net.links) {
        if (down.count(l.id)) continue;
        out[l.from].push_back({l.id, l.to});
        in[l.to].push_back({l.id, l.from});
        if (!inst.net.directed) {
            out[l.to].push_back({l.id, l.from});
            in[l.from].push_back({l.id, l.to});
        }
    }

    for (const auto& dem : inst.demands) {
        const double traffic = (slot < static_cast<int>(dem.v.size())) ? dem.v[slot] : 0.0;
        if (traffic <= 0.0) continue;
        unordered_set<int> must;
        int u = dem.s;
        for (int g = 0; g < N && u != dem.t && out[u].size() == 1; ++g) {
            must.insert(out[u][0].first);
            u = out[u][0].second;
        }
        int v = dem.t;
        for (int g = 0; g < N && v != dem.s && in[v].size() == 1; ++g) {
            must.insert(in[v][0].first);
            v = in[v][0].second;
        }
        for (int lid : must)
            forced[lid] += traffic;
    }
    return forced;
}

// Cached lexicographic objective.
//   unserved  : # demands with positive traffic that could not be routed
//   overflow  : sum of (load - 1) over congested links
//   max_load  : worst link utilization
//   n_crit    : # links in the near-max band (plateau signal: relieving one
//               of several tied hot links is progress before max_load drops)
//   excess    : slot 0 only, sum of (load - tau)^2 over warm links, with
//               tau = kSpreadRel * max_load. Unlike sum_sq it ignores cold
//               links, so a longer route through empty links is not penalized.
//   sum_sq    : sum of squared utilizations (smoothness tie-breaker)
struct Obj {
    int unserved = 0;
    double overflow = 0.0;
    double max_load = 0.0;
    int n_crit = 0;
    double excess = 0.0;
    double sum_sq = 0.0;
};

// Contribution of one link to Obj::excess.
static double excess_term(double load, double tau)
{
    return load > tau ? (load - tau) * (load - tau) : 0.0;
}

// Same unserved / overflow / max_load (the "peak" part of the objective).
static bool same_peak(const Obj& a, const Obj& b)
{
    return a.unserved == b.unserved
        && fabs(a.overflow - b.overflow) <= kEps
        && fabs(a.max_load - b.max_load) <= kEps;
}

// Peak-only comparison (used for predicted candidates, no n_crit).
static bool peak_better(const Obj& a, const Obj& b)
{
    if (a.unserved != b.unserved) return a.unserved < b.unserved;
    if (fabs(a.overflow - b.overflow) > kEps) return a.overflow < b.overflow;
    return a.max_load < b.max_load - kEps;
}

// Lexicographic comparison: smaller is better.
// Priority: unserved > overflow > max_load > n_crit > excess > sum_sq.
static bool better(const Obj& a, const Obj& b)
{
    if (!same_peak(a, b)) return peak_better(a, b);
    if (a.n_crit != b.n_crit) return a.n_crit < b.n_crit;
    if (fabs(a.excess - b.excess) > kEps) return a.excess < b.excess;
    return a.sum_sq < b.sum_sq - kEps;
}

// Records the moments where the accepted solution's max link load dropped.
// Each entry is (global_iteration, new_max_load).
// A global iteration counter is shared across slots so the whole run shows
// up as a single monotone-decreasing series.
// Neighborhood ids used for the history / switch log.
enum Op : int { OP_START = -1, OP_INSERT = 0, OP_CHANGE = 1, OP_DELETE = 2, OP_SHAKE = 3,
                OP_EXTEND = 4, OP_ROBUST = 5, OP_KICK = 6, OP_INSERT2 = 7 };

// Name of a move, as written in the history and statistics.
static const char* op_label(int op)
{
    switch (op) {
        case OP_INSERT: return "insert";
        case OP_CHANGE: return "change";
        case OP_DELETE: return "delete";
        case OP_SHAKE:  return "shake";
        case OP_EXTEND: return "extend";
        case OP_ROBUST: return "robust";
        case OP_KICK:   return "kick";
        case OP_INSERT2: return "insert2";
        default:        return "start";
    }
}

// Moves can be switched off for experiments (ablation studies):
//   SR_DISABLE_OPS=change,shake ./build/solver
// Names: insert insert2 change delete extend shake robust kick, plus "pair" (the
// two-waypoint detours inside the robust phase). insert2 is off by default (in
// the A/B test it left the peak unchanged on average and cost ~3x the time);
// SR_ENABLE_OPS=insert2 turns it on.
static bool op_disabled(const string& name)
{
    auto env_list = [](const char* var) {
        const char* e = getenv(var);
        return e ? "," + string(e) + "," : string();
    };
    static const string off = env_list("SR_DISABLE_OPS");
    static const string on = env_list("SR_ENABLE_OPS");
    auto in = [&](const string& l) { return l.find("," + name + ",") != string::npos; };
    if (name == "insert2" && !in(on)) return true;
    return in(off);
}

struct MaxLoadHistory {
    std::vector<std::pair<long long, double>> entries;
    std::vector<int> entry_ops;                         // move that caused each drop
    std::vector<std::pair<long long, int>> switches;    // (iter, op) when op changes

    // Call after an accepted move. Only records if the new max_load is
    // strictly lower than the last recorded value (or if nothing has been
    // recorded yet).
    void record(long long iter, double max_load, int op)
    {
        if (!std::isfinite(max_load)) return;
        if (entries.empty() || max_load < entries.back().second - 1e-12) {
            entries.emplace_back(iter, max_load);
            entry_ops.push_back(op);
        }
    }

    // Log the active neighborhood; only stores a point when it changes.
    void set_op(long long iter, int op)
    {
        if (switches.empty() || switches.back().second != op)
            switches.emplace_back(iter, op);
    }
};

// Per-slot search state. Owns the candidate waypoints for every demand,
// the aggregated link flow, and the cached objective.
struct SlotSearch {
    const Instance* inst = nullptr;
    const Router* router = nullptr;
    int slot = 0;
    int max_waypoints = 0;     // max #waypoints per demand (segments - 1)
    int n_links = 0;
    vector<double> capacity;   // per link
    vector<double> flow;       // per link: aggregated traffic
    vector<double> forced;     // per link: traffic no routing can avoid
    vector<double> traffic;    // per demand: traffic in this slot
    vector<vector<int>> waypoints;       // current candidate per demand
    vector<vector<int>> prev_waypoints;  // previous slot's waypoints (budget ref)
    vector<Split> links;                 // ECMP link fractions per demand
    vector<char> feasible;               // 1 if route exists
    vector<int> budget_cost;             // #segments changed vs prev, per demand
    int budget_used = 0;                 // sum of budget_cost
    double budget_limit = INF;           // this slot's budget cap
    Obj obj;                             // cached objective
    mt19937 rng;                         // per-slot RNG

    MaxLoadHistory* history = nullptr;   // recorder for max-load drops
    map<string, OpStat>* stats = nullptr; // move statistics (may be null)
    long long iter_index = 0;            // global iteration (set by caller)
    int cur_op = OP_START;               // neighborhood currently running
    int extend_cursor = 0;               // rotating window start for neigh_extend

    // Scratch for rank_by (forward / backward bottleneck DP).
    mutable vector<double> rk_mf, rk_sf, rk_mb, rk_sb;

    // Call right before running a neighborhood (after setting iter_index).
    void begin_op(int op)
    {
        cur_op = op;
        if (history) history->set_op(iter_index, op);
        router->trim_cache();
    }

    // Run one move `f` (returns true if it kept a change) as neighborhood
    // `op`, recording its time and effect on the slot's max load.
    template <class F>
    bool run_op(int op, F f)
    {
        begin_op(op);
        const auto t0 = chrono::steady_clock::now();
        const double before = obj.max_load;
        const bool ok = f();
        if (stats) {
            OpStat& s = (*stats)[op_label(op)];
            ++s.calls;
            s.improved += ok ? 1 : 0;
            s.seconds += chrono::duration<double>(chrono::steady_clock::now() - t0).count();
            s.peak_drop += before - obj.max_load;
        }
        return ok;
    }

    // Initialize for a slot. Warm-starts from `prev` if provided.
    void init(const Instance& in, const Router& rt, int sl,
              const vector<vector<int>>* prev, unsigned seed,
              MaxLoadHistory* hist)
    {
        inst = &in;
        router = &rt;
        slot = sl;
        rng.seed(seed);
        history = hist;
        max_waypoints = max(0, in.scenario.max_segments - 1);
        n_links = rt.max_link_id + 1;

        // Capacity per link (0 means "not a real link").
        capacity.assign(n_links, 0.0);
        for (const auto& l : in.net.links)
            capacity[l.id] = l.capacity;
        flow.assign(n_links, 0.0);

        const int D = static_cast<int>(in.demands.size());
        traffic.resize(D);
        waypoints.assign(D, {});
        prev_waypoints.assign(D, {});
        links.assign(D, {});
        feasible.assign(D, 0);
        budget_cost.assign(D, 0);
        budget_used = 0;
        budget_limit = budget_for_slot(in, sl);

        if (prev)
            prev_waypoints = *prev;

        // Initialize each demand: get its traffic, seed waypoints from the
        // previous slot, and apply the path (accumulating flow + budget).
        for (int d = 0; d < D; ++d) {
            const Demand& dem = in.demands[d];
            traffic[d] = (sl < static_cast<int>(dem.v.size())) ? dem.v[sl] : 0.0;
            if (prev)
                waypoints[d] = (*prev)[d];
            apply_path(d, waypoints[d]);
        }
        forced = forced_flow_at(in, rt.num_nodes, n_links, sl);
        recompute_obj();
    }

    // Load of link i with flow f as the objective sees it: 0 when the link
    // carries only forced traffic, the real load otherwise.
    double eff_load(int i, double f) const
    {
        return f > forced[i] * (1.0 + 1e-9) + 1e-12 ? f / capacity[i] : 0.0;
    }

    // Try to install `w_in` as demand d's waypoints, updating link flow and
    // the segment budget. Returns false (state unchanged) if the budget
    // would be exceeded.
    //
    // Repeated consecutive waypoints, and waypoints equal to the source or
    // target, are dropped first: the checker counts such a "segment" (x, x)
    // in both the budget and max_segments, while routing ignores it.
    bool apply_path(int d, const vector<int>& w_in)
    {
        const Demand& dem = inst->demands[d];
        vector<int> w;
        w.reserve(w_in.size());
        int prev = dem.s;
        for (int x : w_in) {
            if (x == prev) continue;
            w.push_back(x);
            prev = x;
        }
        while (!w.empty() && w.back() == dem.t) w.pop_back();
        const int nc = segment_distance(dem.s, dem.t, prev_waypoints[d], w);
        const int delta = nc - budget_cost[d];
        if (budget_used + delta > budget_limit + 1e-9)
            return false;
        budget_used += delta;
        budget_cost[d] = nc;

        Split new_links;
        const bool ok = router->expand(dem.s, dem.t, w, new_links);
        add_flow(d, -1.0);
        waypoints[d] = w;
        links[d] = std::move(new_links);
        feasible[d] = ok ? 1 : 0;
        add_flow(d, +1.0);
        return true;
    }

    // Recompute the cached objective from current flow + feasibility.
    void recompute_obj()
    {
        obj.unserved = 0;
        obj.overflow = 0.0;
        obj.max_load = 0.0;
        obj.n_crit = 0;
        obj.excess = 0.0;
        obj.sum_sq = 0.0;
        for (size_t d = 0; d < feasible.size(); ++d) {
            if (!feasible[d] && traffic[d] > 0.0)
                ++obj.unserved;
        }
        for (int i = 0; i < n_links; ++i) {
            if (capacity[i] <= 0.0) continue;
            const double load = eff_load(i, flow[i]);
            obj.max_load = max(obj.max_load, load);
            if (load > 1.0)
                obj.overflow += load - 1.0;
            obj.sum_sq += load * load;
        }
        const double thr = crit_threshold(obj.max_load);
        const double tau = spread_tau();
        for (int i = 0; i < n_links; ++i) {
            if (capacity[i] <= 0.0) continue;
            const double load = eff_load(i, flow[i]);
            if (load + 1e-15 >= thr)
                ++obj.n_crit;
            obj.excess += excess_term(load, tau);
        }
    }

    static double crit_threshold(double max_load)
    {
        return max_load - max(1e-9, kCritRel * max(max_load, 1.0));
    }

    // Warm-link threshold for Obj::excess; only slot 0 uses it.
    double spread_tau() const
    {
        return slot == 0 ? kSpreadRel * obj.max_load : INF;
    }

    // Objective of the current state with demand d's flow removed.
    // Precondition: add_flow(d, -1.0) has been applied.
    Obj residual_obj(int d) const
    {
        Obj r;
        r.unserved = obj.unserved - ((!feasible[d] && traffic[d] > 0.0) ? 1 : 0);
        const double tau = spread_tau();
        for (int i = 0; i < n_links; ++i) {
            if (capacity[i] <= 0.0) continue;
            const double load = eff_load(i, flow[i]);
            r.max_load = max(r.max_load, load);
            if (load > 1.0)
                r.overflow += load - 1.0;
            r.excess += excess_term(load, tau);
            r.sum_sq += load * load;
        }
        return r;
    }

    int hottest_link() const
    {
        int best = -1;
        double best_load = -1.0;
        for (int i = 0; i < n_links; ++i) {
            if (capacity[i] <= 0.0) continue;
            const double load = eff_load(i, flow[i]);
            if (load > best_load) {
                best_load = load;
                best = i;
            }
        }
        return best;
    }

    void add_flow(int d, double sign)
    {
        if (!feasible[d] || traffic[d] == 0.0) return;
        for (const auto& [lid, f] : links[d])
            flow[lid] += sign * traffic[d] * f;
    }

    vector<char> critical_mask() const
    {
        vector<char> crit(n_links, 0);
        const double thr = crit_threshold(obj.max_load);
        bool any = false;
        for (int i = 0; i < n_links; ++i) {
            if (capacity[i] <= 0.0) continue;
            if (eff_load(i, flow[i]) + 1e-15 >= thr) {
                crit[i] = 1;
                any = true;
            }
        }
        if (!any) {
            const int hot = hottest_link();
            if (hot >= 0) crit[hot] = 1;
        }
        return crit;
    }

    bool uses_critical(int d, const vector<char>& crit) const
    {
        if (!feasible[d])
            return traffic[d] > 0.0;
        for (const auto& e : links[d]) {
            if (e.first >= 0 && e.first < n_links && crit[e.first])
                return true;
        }
        return false;
    }

    // Demands on the near-max-load cut, largest volume first.
    vector<int> bottleneck_demands() const
    {
        const vector<char> crit = critical_mask();
        vector<int> ds;
        const int D = static_cast<int>(waypoints.size());
        ds.reserve(D);
        for (int d = 0; d < D; ++d) {
            if (traffic[d] <= 0.0 && feasible[d]) continue;
            if (uses_critical(d, crit))
                ds.push_back(d);
        }
        if (ds.empty()) {
            for (int d = 0; d < D; ++d)
                if (traffic[d] > 0.0) ds.push_back(d);
        }
        sort(ds.begin(), ds.end(), [&](int a, int b) {
            if (fabs(traffic[a] - traffic[b]) > 1e-12)
                return traffic[a] > traffic[b];
            return a < b;
        });
        return ds;
    }

    // Raw load of link `link` (0 for non-links).
    double raw_load(int link) const
    {
        return (link < n_links && capacity[link] > 0.0) ? flow[link] / capacity[link] : 0.0;
    }

    // Max and ECMP-weighted sum of link loads over segment a -> b.
    bool sp_load_stats(int a, int b, double& max_load, double& sum_load) const
    {
        max_load = 0.0;
        sum_load = 0.0;
        if (a == b) return true;
        const Split* seg = router->segment(a, b);
        if (!seg) return false;
        for (const auto& [lid, f] : *seg) {
            const double load = raw_load(lid);
            max_load = max(max_load, load);
            sum_load += f * load;
        }
        return true;
    }

    vector<int> route_seq(int d) const
    {
        const Demand& dem = inst->demands[d];
        vector<int> seq;
        seq.reserve(waypoints[d].size() + 2);
        seq.push_back(dem.s);
        for (int w : waypoints[d]) {
            if (w != seq.back())
                seq.push_back(w);
        }
        if (seq.back() != dem.t)
            seq.push_back(dem.t);
        return seq;
    }

    void waypoint_ends(int d, int i, int& pred, int& succ) const
    {
        const Demand& dem = inst->demands[d];
        pred = (i == 0) ? dem.s : waypoints[d][i - 1];
        succ = (i + 1 == static_cast<int>(waypoints[d].size()))
                   ? dem.t : waypoints[d][i + 1];
    }

    int hottest_segment(const vector<int>& seq, double& seg_max, double& seg_sum) const
    {
        int best_i = -1;
        seg_max = -1.0;
        seg_sum = -1.0;
        for (size_t i = 0; i + 1 < seq.size(); ++i) {
            double mx = 0.0, sm = 0.0;
            if (!sp_load_stats(seq[i], seq[i + 1], mx, sm))
                continue;
            if (mx > seg_max + 1e-12 ||
                (fabs(mx - seg_max) <= 1e-12 && sm > seg_sum + 1e-12)) {
                seg_max = mx;
                seg_sum = sm;
                best_i = static_cast<int>(i);
            }
        }
        return best_i;
    }

    bool waypoint_load(int d, int i, double& mx, double& sm) const
    {
        int pred = 0, succ = 0;
        waypoint_ends(d, i, pred, succ);
        const int w = waypoints[d][i];
        double m1 = 0.0, s1 = 0.0, m2 = 0.0, s2 = 0.0;
        if (!sp_load_stats(pred, w, m1, s1) || !sp_load_stats(w, succ, m2, s2))
            return false;
        mx = max(m1, m2);
        sm = s1 + s2;
        return true;
    }

    int hottest_waypoint(int d) const
    {
        int best = -1;
        double best_mx = -1.0, best_sm = -1.0;
        for (int i = 0; i < static_cast<int>(waypoints[d].size()); ++i) {
            double mx = 0.0, sm = 0.0;
            if (!waypoint_load(d, i, mx, sm)) continue;
            if (mx > best_mx + 1e-12 ||
                (fabs(mx - best_mx) <= 1e-12 && sm > best_sm + 1e-12)) {
                best_mx = mx;
                best_sm = sm;
                best = i;
            }
        }
        return best;
    }

    // The k best detour nodes w for a -> w -> b, cheapest first by the max
    // link load `load(link)` over FG(a, w) and FG(w, b), then by the summed
    // load. Both bottlenecks come from one DP each over the shortest-path
    // DAGs (forward from a, backward to b) instead of one walk per w.
    template <class LoadFn>
    vector<int> rank_by(int a, int b, const unordered_set<int>& banned, LoadFn load,
                        int k) const
    {
        const Router& R = *router;
        const int N = R.num_nodes;
        if (a < 0 || b < 0 || a >= N || b >= N) return {};
        rk_mf.assign(N, 0.0);
        rk_sf.assign(N, 0.0);
        rk_mb.assign(N, 0.0);
        rk_sb.assign(N, 0.0);

        // Forward from a: max over FG(a, x) = max over shortest in-arcs p -> x.
        // Nodes in Dijkstra settle order from a.
        const int* from_a = &R.order[R.idx(a, 0)];
        for (int i = 1; i < R.order_len[a]; ++i) {
            const int x = from_a[i];
            const double dx = R.dist[R.idx(a, x)];
            double m = 0.0, s = 0.0;
            int c = 0;
            for (const Arc& e : R.radj[x]) {
                if (R.down[e.link_id]) continue;
                const int p = e.to;
                if (!on_sp(R.dist[R.idx(a, p)], e.metric, dx)) continue;
                const double l = load(e.link_id);
                m = max(m, max(rk_mf[p], l));
                s += rk_sf[p] + l;
                ++c;
            }
            rk_mf[x] = m;
            rk_sf[x] = c ? s / c : 0.0;
        }

        // Backward to b: max over FG(x, b) = max over ECMP next hops x -> y.
        for (int x : R.nodes_to(b)) {
            if (x == b) continue;
            const double dx = R.dist[R.idx(x, b)];
            double m = 0.0, s = 0.0;
            int c = 0;
            for (const Arc& e : R.adj[x]) {
                if (R.down[e.link_id]) continue;
                const int y = e.to;
                if (!on_sp(R.dist[R.idx(y, b)], e.metric, dx)) continue;
                const double l = load(e.link_id);
                m = max(m, max(l, rk_mb[y]));
                s += l + rk_sb[y];
                ++c;
            }
            rk_mb[x] = m;
            rk_sb[x] = c ? s / c : 0.0;
        }

        struct Sc { double mx; double sm; int w; };
        vector<Sc> sc;
        sc.reserve(R.node_ids.size());
        for (int w : R.node_ids) {
            if (banned.count(w) || w == a || w == b) continue;
            if (R.dist[R.idx(a, w)] == INF || R.dist[R.idx(w, b)] == INF) continue;
            sc.push_back({max(rk_mf[w], rk_mb[w]), rk_sf[w] + rk_sb[w], w});
        }
        const auto mid = sc.begin() + min(static_cast<size_t>(max(k, 0)), sc.size());
        partial_sort(sc.begin(), mid, sc.end(), [](const Sc& x, const Sc& y) {
            if (fabs(x.mx - y.mx) > 1e-12) return x.mx < y.mx;
            if (fabs(x.sm - y.sm) > 1e-12) return x.sm < y.sm;
            return x.w < y.w;
        });
        vector<int> out;
        out.reserve(mid - sc.begin());
        for (auto it = sc.begin(); it != mid; ++it) out.push_back(it->w);
        return out;
    }

    // The k best residual detours a -> w -> b, cheapest max-load first.
    vector<int> rank_detour_nodes(int a, int b, const unordered_set<int>& banned, int k) const
    {
        return rank_by(a, b, banned, [&](int lid) { return raw_load(lid); }, k);
    }

    // ---- Delta scoring --------------------------------------------------

    // A fully scored candidate route for one demand.
    struct Cand {
        vector<int> w;
        Obj obj;             // predicted objective (n_crit not predicted)
        int bdelta = 0;      // change in budget_used if applied
        bool hits_crit = false;
    };

    int shortlist_size() const { return slot == 0 ? kShortlistSlot0 : kShortlist; }
    int detour_try() const { return slot == 0 ? kDetourTrySlot0 : kDetourTry; }

    // True if segment a -> b is unreachable or its ECMP split crosses a
    // critical link.
    bool segment_on_cut(int a, int b, const vector<char>& crit) const
    {
        if (a == b) return false;
        const Split* seg = router->segment(a, b);
        if (!seg) return true;
        for (const auto& e : *seg)
            if (e.first < n_links && crit[e.first]) return true;
        return false;
    }

    // Exact objective if demand d were rerouted onto `w`, without applying it.
    // Precondition: d's flow is removed and base == residual_obj(d).
    // Returns false if the move would break the slot's segment budget.
    bool score_candidate(int d, const vector<int>& w, const Obj& base,
                         const vector<char>& crit, Cand& out) const
    {
        const Demand& dem = inst->demands[d];
        out.bdelta = segment_distance(dem.s, dem.t, prev_waypoints[d], w) - budget_cost[d];
        if (budget_used + out.bdelta > budget_limit + 1e-9)
            return false;

        out.w = w;
        out.obj = base;
        out.hits_crit = false;

        Split lk;
        if (!router->expand(dem.s, dem.t, w, lk)) {
            if (traffic[d] > 0.0) ++out.obj.unserved;
            return true;
        }
        if (traffic[d] <= 0.0) return true;

        const double tau = spread_tau();
        for (const auto& [lid, f] : lk) {
            if (lid < 0 || lid >= n_links || capacity[lid] <= 0.0) continue;
            if (crit[lid]) out.hits_crit = true;
            const double r = eff_load(lid, flow[lid]);
            const double n = eff_load(lid, flow[lid] + f * traffic[d]);
            out.obj.max_load = max(out.obj.max_load, n);
            out.obj.overflow += max(0.0, n - 1.0) - max(0.0, r - 1.0);
            out.obj.excess += excess_term(n, tau) - excess_term(r, tau);
            out.obj.sum_sq += n * n - r * r;
        }
        return true;
    }

    // Worth a try_move: lowers the peak, or keeps it while leaving the cut.
    bool promising(const Cand& c) const
    {
        if (peak_better(c.obj, obj)) return true;
        return same_peak(c.obj, obj) && !c.hits_crit;
    }

    // Order: peak first, then "leaves the hot links", then cheaper in
    // segment budget, then smoother. A cheap move that unloads the hottest
    // link beats an expensive one that only helps sum_sq.
    static bool cand_before(const Cand& a, const Cand& b)
    {
        if (!same_peak(a.obj, b.obj)) return peak_better(a.obj, b.obj);
        if (a.hits_crit != b.hits_crit) return !a.hits_crit;
        if (fabs(a.obj.excess - b.obj.excess) > kEps) return a.obj.excess < b.obj.excess;
        if (a.bdelta != b.bdelta) return a.bdelta < b.bdelta;
        return a.obj.sum_sq < b.obj.sum_sq - kEps;
    }

    // Try the best few candidates; keep the first that really improves.
    bool try_ranked(int d, vector<Cand>& cands)
    {
        if (cands.empty()) return false;
        stable_sort(cands.begin(), cands.end(), cand_before);
        const int n = min(detour_try(), static_cast<int>(cands.size()));
        for (int k = 0; k < n; ++k)
            if (try_move(d, cands[k].w))
                return true;
        return false;
    }

    // Shortlist of residual-ranked detour nodes for a -> w -> b.
    vector<int> shortlist(int a, int b, const unordered_set<int>& banned) const
    {
        return rank_detour_nodes(a, b, banned, shortlist_size());
    }

    // ---- Neighborhoods ---------------------------------------------------

    // INSERT: main attack on the bottleneck. For every segment of the route
    // that crosses the critical cut, delta-score a shortlist of detour nodes.
    bool neigh_insert()
    {
        bool any = false;
        vector<char> crit = critical_mask();
        for (int d : bottleneck_demands()) {
            if (static_cast<int>(waypoints[d].size()) >= max_waypoints)
                continue;

            const double peak = obj.max_load;
            const Demand& dem = inst->demands[d];
            add_flow(d, -1.0);
            const Obj base = residual_obj(d);
            const vector<int> seq = route_seq(d);

            unordered_set<int> banned;
            banned.insert(dem.s);
            banned.insert(dem.t);
            for (int w : waypoints[d]) banned.insert(w);

            vector<int> positions;
            for (size_t si = 0; si + 1 < seq.size(); ++si)
                if (segment_on_cut(seq[si], seq[si + 1], crit))
                    positions.push_back(static_cast<int>(si));
            if (positions.empty() && seq.size() >= 2) {
                double smx = 0.0, ssm = 0.0;
                const int si = hottest_segment(seq, smx, ssm);
                if (si >= 0) positions.push_back(si);
            }

            vector<Cand> cands;
            for (int si : positions) {
                for (int node : shortlist(seq[si], seq[si + 1], banned)) {
                    vector<int> w = waypoints[d];
                    w.insert(w.begin() + si, node);
                    Cand c;
                    if (score_candidate(d, w, base, crit, c) && promising(c))
                        cands.push_back(std::move(c));
                }
            }
            add_flow(d, +1.0);

            if (try_ranked(d, cands)) {
                any = true;
                if (obj.max_load < peak - kEps)
                    return true;          // the cut moved: re-target
                crit = critical_mask();
            }
        }
        return any;
    }

    // INSERT2: insert two waypoints at once on a segment that crosses the
    // critical cut. A single waypoint fails when every shortest path to the
    // useful detour node runs through the congested link itself (the node
    // sits behind the link's head, say); a pair (x, y) reaches y from the
    // other side. First waypoints y: the anchors around each critical link
    // u -> v the demand crosses (other in-neighbours of v, other
    // out-neighbours of u) plus the kPairFirst best residual detour nodes.
    // Partners x: the kPairSecond best detour nodes before or after y.
    bool neigh_insert2()
    {
        if (max_waypoints < 2) return false;
        bool any = false;
        vector<char> crit = critical_mask();
        for (int d : bottleneck_demands()) {
            const int nw = static_cast<int>(waypoints[d].size());
            if (nw + 2 > max_waypoints) continue;

            const double peak = obj.max_load;
            const Demand& dem = inst->demands[d];
            const vector<int> seq = route_seq(d);
            unordered_set<int> banned;
            banned.insert(dem.s);
            banned.insert(dem.t);
            for (int w : waypoints[d]) banned.insert(w);

            // Anchors around the critical links on d's route.
            vector<int> anchors;
            for (const auto& e : links[d]) {
                if (e.first >= n_links || !crit[e.first]) continue;
                for (const Arc& a : router->radj[router->link_to[e.first]])
                    if (a.link_id != e.first && !router->down[a.link_id]) anchors.push_back(a.to);
                for (const Arc& a : router->adj[router->link_from[e.first]])
                    if (a.link_id != e.first && !router->down[a.link_id]) anchors.push_back(a.to);
            }

            vector<int> positions;
            for (size_t si = 0; si + 1 < seq.size(); ++si)
                if (segment_on_cut(seq[si], seq[si + 1], crit))
                    positions.push_back(static_cast<int>(si));
            if (positions.empty()) continue;

            add_flow(d, -1.0);
            const Obj base = residual_obj(d);
            vector<Cand> cands;
            for (int si : positions) {
                const int a = seq[si], b = seq[si + 1];
                vector<int> firsts = anchors;
                for (int y : rank_detour_nodes(a, b, banned, kPairFirst)) firsts.push_back(y);
                sort(firsts.begin(), firsts.end());
                firsts.erase(unique(firsts.begin(), firsts.end()), firsts.end());
                for (int y : firsts) {
                    if (banned.count(y) || y == a || y == b) continue;
                    unordered_set<int> banned_y = banned;
                    banned_y.insert(y);
                    // x before y: a -> x -> y -> b
                    for (int x : rank_detour_nodes(a, y, banned_y, kPairSecond)) {
                        vector<int> w = waypoints[d];
                        w.insert(w.begin() + si, {x, y});
                        Cand c;
                        if (score_candidate(d, w, base, crit, c) && promising(c))
                            cands.push_back(std::move(c));
                    }
                    // x after y: a -> y -> x -> b
                    for (int x : rank_detour_nodes(y, b, banned_y, kPairSecond)) {
                        vector<int> w = waypoints[d];
                        w.insert(w.begin() + si, {y, x});
                        Cand c;
                        if (score_candidate(d, w, base, crit, c) && promising(c))
                            cands.push_back(std::move(c));
                    }
                }
            }
            add_flow(d, +1.0);

            if (try_ranked(d, cands)) {
                any = true;
                if (obj.max_load < peak - kEps)
                    return true;          // the cut moved: re-target
                crit = critical_mask();
            }
        }
        return any;
    }

    // Run neighborhood `op` once (insert, insert2, change, delete, extend).
    bool neighborhood(int op)
    {
        switch (op) {
            case OP_INSERT:  return neigh_insert();
            case OP_INSERT2: return neigh_insert2();
            case OP_CHANGE:  return neigh_change();
            case OP_DELETE:  return neigh_delete();
            case OP_EXTEND:  return neigh_extend();
            default:         return false;
        }
    }

    // CHANGE: diversification. Replace any waypoint whose adjacent segments
    // touch the cut with a delta-scored alternative.
    bool neigh_change()
    {
        bool any = false;
        vector<char> crit = critical_mask();
        for (int d : bottleneck_demands()) {
            if (waypoints[d].empty()) continue;

            const double peak = obj.max_load;
            const Demand& dem = inst->demands[d];
            add_flow(d, -1.0);
            const Obj base = residual_obj(d);

            unordered_set<int> banned;
            banned.insert(dem.s);
            banned.insert(dem.t);
            for (int w : waypoints[d]) banned.insert(w);

            vector<Cand> cands;
            const int nw = static_cast<int>(waypoints[d].size());
            for (int i = 0; i < nw; ++i) {
                int pred = 0, succ = 0;
                waypoint_ends(d, i, pred, succ);
                const int wi = waypoints[d][i];
                if (!segment_on_cut(pred, wi, crit) && !segment_on_cut(wi, succ, crit))
                    continue;
                for (int node : shortlist(pred, succ, banned)) {
                    vector<int> w = waypoints[d];
                    w[i] = node;
                    Cand c;
                    if (score_candidate(d, w, base, crit, c) && promising(c))
                        cands.push_back(std::move(c));
                }
            }
            add_flow(d, +1.0);

            if (try_ranked(d, cands)) {
                any = true;
                if (obj.max_load < peak - kEps)
                    return true;
                crit = critical_mask();
            }
        }
        return any;
    }

    // Delete a waypoint only if the residual shortcut is not hotter than
    // the current detour (do not dump the demand back onto the cut).
    bool neigh_delete()
    {
        bool any = false;
        const int D = static_cast<int>(waypoints.size());
        vector<int> ds = bottleneck_demands();
        vector<char> seen(D, 0);
        for (int d : ds) seen[d] = 1;
        for (int d = 0; d < D; ++d) {
            if (!seen[d] && !waypoints[d].empty() && traffic[d] > 0.0)
                ds.push_back(d);
        }

        for (int d : ds) {
            if (waypoints[d].empty()) continue;

            add_flow(d, -1.0);
            int best_i = -1;
            double best_gain = -INF;
            for (int i = 0; i < static_cast<int>(waypoints[d].size()); ++i) {
                int pred = 0, succ = 0;
                waypoint_ends(d, i, pred, succ);
                double m1 = 0.0, s1 = 0.0, m2 = 0.0, s2 = 0.0, ms = 0.0, ss = 0.0;
                if (!sp_load_stats(pred, waypoints[d][i], m1, s1) ||
                    !sp_load_stats(waypoints[d][i], succ, m2, s2) ||
                    !sp_load_stats(pred, succ, ms, ss))
                    continue;
                const double detour = max(m1, m2);
                if (ms > detour + 1e-12) continue; // shortcut sits on a hotter residual
                const double gain = detour - ms;
                if (gain > best_gain + 1e-12) {
                    best_gain = gain;
                    best_i = i;
                }
            }
            add_flow(d, +1.0);
            if (best_i < 0) continue;

            vector<int> cand = waypoints[d];
            cand.erase(cand.begin() + best_i);
            if (try_move(d, cand))
                any = true;
        }
        return any;
    }

    // Worth a try_move for EXTEND: lowers the peak, or keeps it while
    // moving load off warm links.
    bool spreads(const Cand& c) const
    {
        if (peak_better(c.obj, obj)) return true;
        return same_peak(c.obj, obj) && c.obj.excess < obj.excess - kEps;
    }

    // EXTEND (slot 0 only): lengthen routes through the cold part of the
    // network. For demands crossing warm links (load >= kSpreadRel * max_load),
    // insert one more waypoint on a warm segment if that keeps the peak and
    // lowers the warm-link excess. Repeated calls grow a route one waypoint at
    // a time up to max_segments; neigh_delete trims waypoints that stop paying.
    bool neigh_extend()
    {
        if (slot != 0 || max_waypoints == 0) return false;

        const int D = static_cast<int>(waypoints.size());
        double tau = spread_tau();
        vector<int> ds;
        for (int d = 0; d < D; ++d) {
            if (traffic[d] <= 0.0 || !feasible[d]) continue;
            if (static_cast<int>(waypoints[d].size()) >= max_waypoints) continue;
            for (const auto& e : links[d]) {
                const int lid = e.first;
                if (lid < n_links && capacity[lid] > 0.0
                    && eff_load(lid, flow[lid]) >= tau) {
                    ds.push_back(d);
                    break;
                }
            }
        }
        if (ds.empty()) return false;
        sort(ds.begin(), ds.end(), [&](int a, int b) {
            if (fabs(traffic[a] - traffic[b]) > 1e-12)
                return traffic[a] > traffic[b];
            return a < b;
        });

        // Rotate over the list so successive calls cover different demands.
        const int n_ds = static_cast<int>(ds.size());
        const int n = min(kExtendDemands, n_ds);
        const int start = extend_cursor % n_ds;
        extend_cursor = (start + n) % n_ds;

        bool any = false;
        for (int k = 0; k < n; ++k) {
            const int d = ds[(start + k) % n_ds];
            if (static_cast<int>(waypoints[d].size()) >= max_waypoints) continue;

            tau = spread_tau();
            const vector<char> crit = critical_mask();
            const Demand& dem = inst->demands[d];
            const vector<int> seq = route_seq(d);

            // Warm segments, measured with d's own flow still on the network.
            vector<int> positions;
            for (size_t si = 0; si + 1 < seq.size(); ++si) {
                double mx = 0.0, sm = 0.0;
                if (sp_load_stats(seq[si], seq[si + 1], mx, sm) && mx >= tau)
                    positions.push_back(static_cast<int>(si));
            }
            if (positions.empty()) continue;

            add_flow(d, -1.0);
            const Obj base = residual_obj(d);

            unordered_set<int> banned;
            banned.insert(dem.s);
            banned.insert(dem.t);
            for (int w : waypoints[d]) banned.insert(w);

            vector<Cand> cands;
            for (int si : positions) {
                for (int node : shortlist(seq[si], seq[si + 1], banned)) {
                    vector<int> w = waypoints[d];
                    w.insert(w.begin() + si, node);
                    Cand c;
                    if (score_candidate(d, w, base, crit, c) && spreads(c))
                        cands.push_back(std::move(c));
                }
            }
            add_flow(d, +1.0);

            if (try_ranked(d, cands))
                any = true;
        }
        return any;
    }

    void shake(int strength)
    {
        const vector<int> ds = bottleneck_demands();
        if (ds.empty() || strength <= 0) return;
        uniform_int_distribution<int> pick(0, static_cast<int>(ds.size()) - 1);

        for (int s = 0; s < strength; ++s) {
            const int d = ds[pick(rng)];
            const Demand& dem = inst->demands[d];
            vector<int> cand = waypoints[d];

            int mode = static_cast<int>(rng() % 3);
            if (mode == 0 && static_cast<int>(cand.size()) >= max_waypoints) mode = 1;
            if ((mode == 1 || mode == 2) && cand.empty()) mode = 0;
            if (mode == 0 && static_cast<int>(cand.size()) >= max_waypoints)
                continue;

            add_flow(d, -1.0);
            bool have = false;
            if (mode == 0) {
                const vector<int> seq = route_seq(d);
                if (seq.size() >= 2) {
                    double smx = 0.0, ssm = 0.0;
                    int si = hottest_segment(seq, smx, ssm);
                    if (si < 0) si = 0;
                    unordered_set<int> banned;
                    banned.insert(dem.s);
                    banned.insert(dem.t);
                    for (int w : cand) banned.insert(w);
                    const vector<int> ranked =
                        rank_detour_nodes(seq[si], seq[si + 1], banned, kDetourRank);
                    add_flow(d, +1.0);
                    if (!ranked.empty()) {
                        const int n = min(kDetourRank, static_cast<int>(ranked.size()));
                        uniform_int_distribution<int> rw(0, n - 1);
                        cand.insert(cand.begin() + si, ranked[rw(rng)]);
                        have = true;
                    }
                } else {
                    add_flow(d, +1.0);
                }
            } else if (mode == 1) {
                int i = hottest_waypoint(d);
                if (i < 0)
                    i = static_cast<int>(rng() % cand.size());
                int pred = 0, succ = 0;
                waypoint_ends(d, i, pred, succ);
                unordered_set<int> banned;
                banned.insert(dem.s);
                banned.insert(dem.t);
                for (int w : cand) banned.insert(w);
                const vector<int> ranked = rank_detour_nodes(pred, succ, banned, kDetourRank);
                add_flow(d, +1.0);
                if (!ranked.empty()) {
                    const int n = min(kDetourRank, static_cast<int>(ranked.size()));
                    uniform_int_distribution<int> rw(0, n - 1);
                    cand[i] = ranked[rw(rng)];
                    have = true;
                }
            } else {
                int i = hottest_waypoint(d);
                add_flow(d, +1.0);
                if (i < 0)
                    i = static_cast<int>(rng() % cand.size());
                cand.erase(cand.begin() + i);
                have = true;
            }

            if (have)
                commit_move(d, cand, Accept::Shake);
        }
    }

    enum class Accept { Improve, Shake };

    bool try_move(int d, const vector<int>& cand)
    {
        return commit_move(d, cand, Accept::Improve);
    }

    // Acceptance for descent moves.
    //  - strict lexicographic improvement, except that in later slots a move
    //    that only improves sum_sq may not spend extra segment budget;
    //  - budget recovery: same objective, fewer changed segments.
    bool improves(const Obj& n, const Obj& o, int bdelta) const
    {
        if (better(n, o)) {
            if (slot > 0 && bdelta > 0 && same_peak(n, o) && n.n_crit == o.n_crit)
                return false;
            return true;
        }
        return bdelta < 0 && same_peak(n, o) && n.n_crit <= o.n_crit
            && n.excess <= o.excess + kEps && n.sum_sq <= o.sum_sq + kEps;
    }

    // Acceptance for shake moves: may worsen the peak slightly, never
    // un-serve a demand.
    static bool shake_ok(const Obj& n, const Obj& o)
    {
        return n.unserved <= o.unserved
            && n.overflow <= o.overflow + kShakeSlack
            && n.max_load <= o.max_load * (1.0 + kShakeSlack) + kEps;
    }

    // Everything a single-demand move can change, for undo.
    struct DemandState {
        vector<int> w;
        Split links;
        char feasible = 0;
        int bcost = 0;
        int bused = 0;
        Obj obj;
    };

    DemandState keep(int d) const
    {
        return {waypoints[d], links[d], feasible[d], budget_cost[d], budget_used, obj};
    }

    void undo(int d, const DemandState& s)
    {
        if (feasible[d]) {
            for (const auto& [lid, f] : links[d])
                flow[lid] -= traffic[d] * f;
        }
        waypoints[d] = s.w;
        links[d] = s.links;
        feasible[d] = s.feasible;
        budget_cost[d] = s.bcost;
        budget_used = s.bused;
        if (s.feasible) {
            for (const auto& [lid, f] : links[d])
                flow[lid] += traffic[d] * f;
        }
        obj = s.obj;
    }

    bool commit_move(int d, const vector<int>& cand, Accept mode)
    {
        const DemandState saved = keep(d);
        const Obj old_obj = obj;

        if (!apply_path(d, cand))
            return false;

        recompute_obj();
        const bool keep_it = (mode == Accept::Shake)
                                 ? shake_ok(obj, old_obj)
                                 : improves(obj, old_obj, budget_used - saved.bused);
        if (keep_it) {
            if (history)
                history->record(iter_index, obj.max_load, cur_op); // drops only
            return true;
        }
        undo(d, saved);
        return false;
    }

    // ---- Best-solution snapshot (VNS: shake from the best basin) ----------
    struct Snapshot {
        vector<vector<int>> waypoints;
        vector<Split> links;
        vector<char> feasible;
        vector<int> budget_cost;
        int budget_used = 0;
        vector<double> flow;
        Obj obj;
    };

    void save(Snapshot& s) const
    {
        s.waypoints = waypoints;
        s.links = links;
        s.feasible = feasible;
        s.budget_cost = budget_cost;
        s.budget_used = budget_used;
        s.flow = flow;
        s.obj = obj;
    }

    void restore(const Snapshot& s)
    {
        waypoints = s.waypoints;
        links = s.links;
        feasible = s.feasible;
        budget_cost = s.budget_cost;
        budget_used = s.budget_used;
        flow = s.flow;
        obj = s.obj;
    }
};

// ---------------------------------------------------------------------------
// Robust slot 0
// ---------------------------------------------------------------------------
//
// Budgets are tiny next to the number of demands, so the slot-0 waypoints are
// in practice the waypoints of every slot. FutureSlots keeps the link flow of
// each later slot t as if the slot-0 waypoints stayed in place, under that
// slot's traffic and failed links.
//
// Failures only remove links, so a segment keeps its slot-0 ECMP split at t
// unless its slot-0 forwarding graph uses a link that is down at t. Only those
// segments are re-split, with distances to the segment end repaired from the
// slot-0 ones (only nodes whose shortest paths all crossed a failed link are
// re-solved).

struct FutureSlots {
    struct FailedArc { int tail; int head; double metric; };
    // Distances to one node at slot t that differ from slot 0, sorted by node.
    struct DistDiff { vector<pair<int, double>> changed; };
    struct SegT { bool ok = false; Split s; };

    const Instance* inst = nullptr;
    const Router* r0 = nullptr;     // slot-0 router
    int T = 0;
    int N = 0;
    int n_links = 0;
    vector<double> capacity;
    vector<uint64_t> link_mask;          // bit t set if the link is down at t (t < 64)
    vector<vector<char>> down;           // [t][link]
    vector<vector<FailedArc>> failed;    // [t] arcs down at t but up at slot 0
    vector<char> full;                   // [t] a slot-0 failure is repaired: re-solve fully
    vector<vector<double>> traffic;      // [t][demand]
    vector<vector<double>> forced;       // [t][link]
    vector<vector<double>> flow;         // [t][link]
    vector<int> unserved;                // [t]

    unordered_map<uint64_t, DistDiff> rows;     // (t, end) -> distance diff
    unordered_map<uint64_t, SegT> segs;         // (t, a, b) -> split at t
    unordered_map<uint64_t, uint64_t> masks;    // (a, b) -> slots that touch FG0(a, b)
    size_t cached = 0;
    EcmpScratch sc;
    vector<int> mark, inc, fin;                 // stamps for row repair
    vector<double> key;
    int stamp = 0;

    void init(const Instance& in, const Router& r, const SlotSearch& st)
    {
        inst = &in;
        r0 = &r;
        T = in.num_time_slots;
        N = r.num_nodes;
        n_links = r.max_link_id + 1;
        capacity = st.capacity;
        link_mask.assign(n_links, 0);
        down.assign(T, vector<char>(n_links, 0));
        failed.assign(T, {});
        full.assign(T, 0);
        for (const auto& iv : in.scenario.interventions) {
            if (iv.t < 1 || iv.t >= T) continue;
            for (int lid : iv.links)
                if (lid >= 0 && lid < n_links) down[iv.t][lid] = 1;
        }
        for (int t = 1; t < T; ++t) {
            for (const auto& l : in.net.links) {
                if (down[t][l.id]) {
                    if (t < 64) link_mask[l.id] |= uint64_t(1) << t;
                    if (r.down[l.id]) continue;
                    failed[t].push_back({l.from, l.to, l.metric});
                    if (!in.net.directed)
                        failed[t].push_back({l.to, l.from, l.metric});
                } else if (r.down[l.id]) {
                    full[t] = 1;
                }
            }
        }

        const int D = static_cast<int>(in.demands.size());
        traffic.assign(T, vector<double>(D, 0.0));
        forced.assign(T, {});
        flow.assign(T, vector<double>(n_links, 0.0));
        unserved.assign(T, 0);
        for (int t = 1; t < T; ++t) {
            for (int d = 0; d < D; ++d) {
                const auto& v = in.demands[d].v;
                traffic[t][d] = (t < static_cast<int>(v.size())) ? v[t] : 0.0;
            }
            forced[t] = forced_flow_at(in, N, n_links, t);
        }
        mark.assign(N, 0);
        inc.assign(N, 0);
        fin.assign(N, 0);
        key.assign(N, INF);
        for (int d = 0; d < D; ++d)
            apply(d, st.waypoints[d], +1.0);
    }

    double eff_load(int t, int i, double f) const
    {
        return f > forced[t][i] * (1.0 + 1e-9) + 1e-12 ? f / capacity[i] : 0.0;
    }

    void trim()
    {
        if (cached > kFutureCacheCap) {
            rows.clear();
            segs.clear();
            masks.clear();
            cached = 0;
        }
    }

    double dist_t(const DistDiff& dd, int x, int b) const
    {
        auto it = lower_bound(dd.changed.begin(), dd.changed.end(),
                              make_pair(x, -INF));
        if (it != dd.changed.end() && it->first == x) return it->second;
        return r0->dist[r0->idx(x, b)];
    }

    // Distances to b at slot t, as a diff against slot 0.
    const DistDiff& row(int t, int b)
    {
        const uint64_t k = (static_cast<uint64_t>(t) << 32) | static_cast<uint32_t>(b);
        auto it = rows.find(k);
        if (it != rows.end()) return it->second;
        DistDiff dd;
        repair_row(t, b, dd);
        cached += dd.changed.size() + 1;
        return rows.emplace(k, std::move(dd)).first->second;
    }

    void repair_row(int t, int b, DistDiff& dd)
    {
        auto D0 = [&](int x) { return r0->dist[r0->idx(x, b)]; };
        using P = pair<double, int>;

        if (full[t]) {
            // Plain reverse Dijkstra on the slot-t topology.
            vector<double> d(N, INF);
            priority_queue<P, vector<P>, greater<P>> pq;
            d[b] = 0.0;
            pq.push({0.0, b});
            while (!pq.empty()) {
                auto [dx, x] = pq.top();
                pq.pop();
                if (dx > d[x]) continue;
                for (const Arc& e : r0->radj[x]) {
                    if (down[t][e.link_id]) continue;
                    if (dx + e.metric < d[e.to]) {
                        d[e.to] = dx + e.metric;
                        pq.push({d[e.to], e.to});
                    }
                }
            }
            for (int x : r0->node_ids) dd.changed.push_back({x, d[x]});
            sort(dd.changed.begin(), dd.changed.end());
            return;
        }

        if (++stamp == numeric_limits<int>::max()) {
            fill(mark.begin(), mark.end(), 0);
            fill(inc.begin(), inc.end(), 0);
            fill(fin.begin(), fin.end(), 0);
            stamp = 1;
        }
        const int s = stamp;

        // Nodes whose shortest paths to b may cross a failed arc: tails of
        // failed shortest-path arcs and their shortest-path ancestors.
        vector<int> cand;
        for (const FailedArc& f : failed[t]) {
            if (mark[f.tail] != s && on_sp(D0(f.head), f.metric, D0(f.tail))) {
                mark[f.tail] = s;
                cand.push_back(f.tail);
            }
        }
        if (cand.empty()) return;
        for (size_t i = 0; i < cand.size(); ++i) {
            const int x = cand[i];
            const double dx = D0(x);
            for (const Arc& e : r0->radj[x]) {
                const int z = e.to;
                if (r0->down[e.link_id] || mark[z] == s) continue;
                if (!on_sp(dx, e.metric, D0(z))) continue;
                mark[z] = s;
                cand.push_back(z);
            }
        }

        // A candidate keeps its distance if one usable next hop does.
        sort(cand.begin(), cand.end(), [&](int a, int c) { return D0(a) < D0(c); });
        vector<int> raised;
        for (int x : cand) {
            bool kept = false;
            for (const Arc& e : r0->adj[x]) {
                if (down[t][e.link_id] || inc[e.to] == s) continue;
                if (on_sp(D0(e.to), e.metric, D0(x))) { kept = true; break; }
            }
            if (!kept) {
                inc[x] = s;
                raised.push_back(x);
            }
        }
        if (raised.empty()) return;

        // Re-solve the raised nodes from their unchanged neighbours.
        priority_queue<P, vector<P>, greater<P>> pq;
        for (int x : raised) {
            key[x] = INF;
            for (const Arc& e : r0->adj[x]) {
                if (down[t][e.link_id] || inc[e.to] == s) continue;
                const double dy = D0(e.to);
                if (dy < INF) key[x] = min(key[x], dy + e.metric);
            }
            if (key[x] < INF) pq.push({key[x], x});
        }
        while (!pq.empty()) {
            auto [kx, x] = pq.top();
            pq.pop();
            if (fin[x] == s || kx > key[x]) continue;
            fin[x] = s;
            for (const Arc& e : r0->radj[x]) {
                const int z = e.to;
                if (down[t][e.link_id] || inc[z] != s || fin[z] == s) continue;
                if (kx + e.metric < key[z]) {
                    key[z] = kx + e.metric;
                    pq.push({key[z], z});
                }
            }
        }
        for (int x : raised)
            dd.changed.push_back({x, fin[x] == s ? key[x] : INF});
        sort(dd.changed.begin(), dd.changed.end());
    }

    // True if FG0(a, b) (split s0) uses a link that is down at t.
    bool affected(int t, int a, int b, const Split& s0)
    {
        if (t < 64) {
            const uint64_t k = pack_seg(a, b);
            auto it = masks.find(k);
            uint64_t m = 0;
            if (it != masks.end()) {
                m = it->second;
            } else {
                for (const auto& e : s0) m |= link_mask[e.first];
                masks.emplace(k, m);
                ++cached;
            }
            return (m >> t) & 1;
        }
        for (const auto& e : s0)
            if (down[t][e.first]) return true;
        return false;
    }

    // ECMP split of segment a -> b at slot t; nullptr if unreachable.
    const Split* seg(int t, int a, int b)
    {
        static const Split kNoLinks;
        if (a == b) return &kNoLinks;
        if (!full[t]) {
            const Split* s0 = r0->segment(a, b);
            if (!s0) return nullptr;     // failures never reconnect a pair
            if (!affected(t, a, b, *s0)) return s0;
        }
        const uint64_t k = (static_cast<uint64_t>(t) << 48)
                         | (static_cast<uint64_t>(a) << 24) | static_cast<uint64_t>(b);
        auto it = segs.find(k);
        if (it == segs.end()) {
            SegT st;
            const DistDiff& dd = row(t, b);
            if (dist_t(dd, a, b) < INF) {
                st.ok = true;
                ecmp_split(a, b, r0->adj, down[t],
                           [&](int x) { return dist_t(dd, x, b); }, sc, st.s);
            }
            cached += st.s.size() + 1;
            it = segs.emplace(k, std::move(st)).first;
        }
        return it->second.ok ? &it->second.s : nullptr;
    }

    // Link fractions of demand d routed on waypoints w at slot t.
    bool route(int t, int d, const vector<int>& w, Split& out)
    {
        const Demand& dem = inst->demands[d];
        out.clear();
        int prev = dem.s;
        auto add = [&](int a, int b) {
            if (a == b) return true;
            const Split* s = seg(t, a, b);
            if (!s) return false;
            out.insert(out.end(), s->begin(), s->end());
            return true;
        };
        for (int x : w) {
            if (x == prev) continue;
            if (!add(prev, x)) { out.clear(); return false; }
            prev = x;
        }
        if (!add(prev, dem.t)) { out.clear(); return false; }
        merge_split(out);
        return true;
    }

    // Add (sign = +1) or remove (-1) demand d on waypoints w in every later slot.
    void apply(int d, const vector<int>& w, double sign)
    {
        Split lk;
        for (int t = 1; t < T; ++t) {
            const double tr = traffic[t][d];
            if (tr <= 0.0) continue;
            if (route(t, d, w, lk)) {
                for (const auto& [lid, f] : lk)
                    flow[t][lid] += sign * tr * f;
            } else {
                unserved[t] += (sign > 0.0) ? 1 : -1;
            }
        }
    }
};

// Objective over all slots: slot 0 (current state) plus every later slot
// (FutureSlots). Same lexicographic order as Obj, without sum_sq.
struct GObj {
    int unserved = 0;
    double overflow = 0.0;
    double peak = 0.0;
    int n_crit = 0;          // (slot, link) pairs in the near-peak band
    double excess = 0.0;     // sum of (load - tau)^2 over warm (slot, link) pairs
};

static bool g_same_peak(const GObj& a, const GObj& b)
{
    return a.unserved == b.unserved
        && fabs(a.overflow - b.overflow) <= kEps
        && fabs(a.peak - b.peak) <= kEps;
}

static bool g_peak_better(const GObj& a, const GObj& b)
{
    if (a.unserved != b.unserved) return a.unserved < b.unserved;
    if (fabs(a.overflow - b.overflow) > kEps) return a.overflow < b.overflow;
    return a.peak < b.peak - kEps;
}

static bool g_better(const GObj& a, const GObj& b)
{
    if (!g_same_peak(a, b)) return g_peak_better(a, b);
    if (a.n_crit != b.n_crit) return a.n_crit < b.n_crit;
    return a.excess < b.excess - kEps;
}

// ROBUST (slot 0 only): like EXTEND, but every slot counts. Each link gets a
// profile = its worst load over all slots; demands crossing warm profile
// links get a waypoint inserted, or one replaced, through detour nodes whose
// forwarding graphs have the lowest profile ("fill the low-load nodes").
// Candidates are scored exactly over all slots; a move is kept only if it
// lowers the all-slot objective (peak, then near-peak count, then excess).
struct RobustSearch {
    SlotSearch& st;
    FutureSlots& fut;
    GObj cur;
    double tau = 0.0;
    vector<double> profile;     // per link: max effective load over all slots
    double thr = 0.0;           // near-peak threshold of cur.peak
    vector<char> hot_slot;      // per slot: its peak is in the near-peak band
    vector<int> link_from, link_to;  // link id -> endpoints
    int hot_cursor = 0;
    int warm_cursor = 0;
    int warm_size = 0;          // size of the last warm list (for stopping)
    vector<int> last_hot, last_warm;  // demand lists of the last step (for kicks)

    RobustSearch(SlotSearch& s, FutureSlots& f) : st(s), fut(f)
    {
        link_from.assign(st.n_links, -1);
        link_to.assign(st.n_links, -1);
        for (const auto& l : st.inst->net.links) {
            link_from[l.id] = l.from;
            link_to[l.id] = l.to;
        }
    }

    double load_at(int t, int i) const
    {
        return t == 0 ? st.eff_load(i, st.flow[i]) : fut.eff_load(t, i, fut.flow[t][i]);
    }

    int base_unserved() const
    {
        int u = st.obj.unserved;
        for (int t = 1; t < fut.T; ++t) u += fut.unserved[t];
        return u;
    }

    // Full all-slot objective with warm threshold tau_ref (< 0: own peak).
    GObj evaluate(double tau_ref) const
    {
        GObj g;
        g.unserved = base_unserved();
        for (int t = 0; t < fut.T; ++t) {
            for (int i = 0; i < st.n_links; ++i) {
                if (st.capacity[i] <= 0.0) continue;
                const double load = load_at(t, i);
                g.peak = max(g.peak, load);
                if (load > 1.0) g.overflow += load - 1.0;
            }
        }
        const double tr = tau_ref < 0.0 ? kSpreadRel * g.peak : tau_ref;
        const double thr = SlotSearch::crit_threshold(g.peak);
        for (int t = 0; t < fut.T; ++t) {
            for (int i = 0; i < st.n_links; ++i) {
                if (st.capacity[i] <= 0.0) continue;
                const double load = load_at(t, i);
                if (load + 1e-15 >= thr) ++g.n_crit;
                g.excess += excess_term(load, tr);
            }
        }
        return g;
    }

    void refresh()
    {
        cur = evaluate(-1.0);
        tau = kSpreadRel * cur.peak;
    }

    // All-slot objective with demand d removed everywhere (no n_crit).
    GObj residual(int d) const
    {
        GObj g;
        g.unserved = base_unserved()
                   - ((!st.feasible[d] && st.traffic[d] > 0.0) ? 1 : 0);
        for (int t = 0; t < fut.T; ++t) {
            for (int i = 0; i < st.n_links; ++i) {
                if (st.capacity[i] <= 0.0) continue;
                const double load = load_at(t, i);
                g.peak = max(g.peak, load);
                if (load > 1.0) g.overflow += load - 1.0;
                g.excess += excess_term(load, tau);
            }
        }
        return g;
    }

    // Exact all-slot objective if d were routed on w (d already removed).
    GObj score(int d, const vector<int>& w, const GObj& base)
    {
        GObj g = base;
        Split lk;
        auto add = [&](int t, double tr, const vector<double>& fl) {
            for (const auto& [lid, f] : lk) {
                if (lid < 0 || lid >= st.n_links || st.capacity[lid] <= 0.0) continue;
                const double r = t == 0 ? st.eff_load(lid, fl[lid])
                                        : fut.eff_load(t, lid, fl[lid]);
                const double n = t == 0 ? st.eff_load(lid, fl[lid] + f * tr)
                                        : fut.eff_load(t, lid, fl[lid] + f * tr);
                g.peak = max(g.peak, n);
                g.overflow += max(0.0, n - 1.0) - max(0.0, r - 1.0);
                g.excess += excess_term(n, tau) - excess_term(r, tau);
            }
        };
        const Demand& dem = st.inst->demands[d];
        if (st.traffic[d] > 0.0) {
            if (st.router->expand(dem.s, dem.t, w, lk)) add(0, st.traffic[d], st.flow);
            else ++g.unserved;
        }
        for (int t = 1; t < fut.T; ++t) {
            const double tr = fut.traffic[t][d];
            if (tr <= 0.0) continue;
            if (fut.route(t, d, w, lk)) add(t, tr, fut.flow[t]);
            else ++g.unserved;
        }
        return g;
    }

    bool worth(const GObj& c) const
    {
        if (g_peak_better(c, cur)) return true;
        return g_same_peak(c, cur) && c.excess < cur.excess - kEps;
    }

    // Install w for d in slot 0 and every later slot; keep it only if the
    // all-slot objective improves.
    bool commit(int d, const vector<int>& w)
    {
        const SlotSearch::DemandState saved = st.keep(d);
        const vector<int> old_w = st.waypoints[d];
        if (!st.apply_path(d, w))
            return false;
        st.recompute_obj();
        fut.apply(d, old_w, -1.0);
        fut.apply(d, w, +1.0);
        const GObj n = evaluate(tau);
        if (g_better(n, cur)) {
            refresh();
            if (st.history)
                st.history->record(st.iter_index, st.obj.max_load, OP_ROBUST);
            return true;
        }
        fut.apply(d, w, -1.0);
        fut.apply(d, old_w, +1.0);
        st.undo(d, saved);
        return false;
    }

    // Worst load on segment a -> b over all slots, each slot with its own
    // split (a failure can re-route the segment onto other links).
    double seg_heat(int a, int b)
    {
        if (a == b) return 0.0;
        double m = 0.0;
        if (const Split* s = st.router->segment(a, b))
            for (const auto& e : *s) m = max(m, load_at(0, e.first));
        for (int t = 1; t < fut.T; ++t)
            if (const Split* s = fut.seg(t, a, b))
                for (const auto& e : *s) m = max(m, load_at(t, e.first));
        return m;
    }

    vector<int> shortlist(int a, int b, const unordered_set<int>& banned, int width) const
    {
        return st.rank_by(a, b, banned, [&](int lid) { return profile[lid]; }, width);
    }

    // Nodes that lead around a near-peak link u -> v crossed by d at a hot
    // slot: tails of v's other usable in-links and heads of u's other usable
    // out-links. Reaching such a node from the far side often needs a second
    // waypoint, so they are tried alone and in pairs (see improve_demand).
    vector<int> detour_anchors(int d)
    {
        vector<int> out;
        Split lk;
        for (int t = 0; t < fut.T; ++t) {
            if (!hot_slot[t]) continue;
            if (t == 0) {
                if (st.traffic[d] <= 0.0) continue;
                lk = st.links[d];
            } else if (fut.traffic[t][d] <= 0.0 || !fut.route(t, d, st.waypoints[d], lk)) {
                continue;
            }
            const vector<char>& down = t == 0 ? st.router->down : fut.down[t];
            for (const auto& e : lk) {
                if (load_at(t, e.first) + 1e-15 < thr) continue;
                for (const Arc& a : st.router->radj[link_to[e.first]])
                    if (a.link_id != e.first && !down[a.link_id]) out.push_back(a.to);
                for (const Arc& a : st.router->adj[link_from[e.first]])
                    if (a.link_id != e.first && !down[a.link_id]) out.push_back(a.to);
            }
        }
        sort(out.begin(), out.end());
        out.erase(unique(out.begin(), out.end()), out.end());
        return out;
    }

    struct RCand { vector<int> w; GObj g; };

    // Candidate routes for d, scored exactly over all slots; those `keep(g)`
    // accepts are returned. Moves: insert or replace one waypoint through
    // low-profile nodes (`width` per position: the profile ranking uses
    // slot-0 paths, so for demands on a later slot's peak, re-routed by a
    // failure, the useful nodes can sit far down the list), and for hot
    // demands two-waypoint detours around the peak links.
    template <class Keep>
    vector<RCand> candidates(int d, int width, bool hot, Keep keep)
    {
        const Demand& dem = st.inst->demands[d];
        const vector<int> seq = st.route_seq(d);
        const vector<int> cur_w = st.waypoints[d];

        vector<int> positions;
        if (static_cast<int>(cur_w.size()) < st.max_waypoints)
            for (size_t si = 0; si + 1 < seq.size(); ++si)
                if (seg_heat(seq[si], seq[si + 1]) >= tau)
                    positions.push_back(static_cast<int>(si));
        vector<int> changes;
        for (int i = 0; i < static_cast<int>(cur_w.size()); ++i) {
            int pred = 0, succ = 0;
            st.waypoint_ends(d, i, pred, succ);
            if (seg_heat(pred, cur_w[i]) >= tau || seg_heat(cur_w[i], succ) >= tau)
                changes.push_back(i);
        }
        vector<RCand> cands;
        if (positions.empty() && changes.empty() && !hot) return cands;

        unordered_set<int> banned;
        banned.insert(dem.s);
        banned.insert(dem.t);
        for (int w : cur_w) banned.insert(w);

        static const bool pairs_on = !op_disabled("pair");
        const vector<int> anchors = hot && pairs_on ? detour_anchors(d) : vector<int>{};
        st.add_flow(d, -1.0);
        fut.apply(d, cur_w, -1.0);
        const GObj base = residual(d);

        for (int si : positions) {
            for (int node : shortlist(seq[si], seq[si + 1], banned, width)) {
                vector<int> w = cur_w;
                w.insert(w.begin() + si, node);
                GObj g = score(d, w, base);
                if (keep(g)) cands.push_back({std::move(w), g});
            }
        }
        for (int i : changes) {
            int pred = 0, succ = 0;
            st.waypoint_ends(d, i, pred, succ);
            for (int node : shortlist(pred, succ, banned, width)) {
                vector<int> w = cur_w;
                w[i] = node;
                GObj g = score(d, w, base);
                if (keep(g)) cands.push_back({std::move(w), g});
            }
        }
        // Two-waypoint detours around the peak (hot demands only).
        const int nw = static_cast<int>(cur_w.size());
        for (int y : (hot && nw < st.max_waypoints) ? anchors : vector<int>{}) {
            if (banned.count(y)) continue;
            unordered_set<int> banned_y = banned;
            banned_y.insert(y);
            for (int pos = 0; pos <= nw; ++pos) {
                vector<int> w1 = cur_w;
                w1.insert(w1.begin() + pos, y);
                GObj g = score(d, w1, base);
                if (keep(g)) cands.push_back({w1, g});
                if (nw + 2 > st.max_waypoints) continue;
                const int prev = pos == 0 ? dem.s : cur_w[pos - 1];
                const int next = pos == nw ? dem.t : cur_w[pos];
                for (int x : shortlist(prev, y, banned_y, kRobustPair)) {
                    vector<int> w2 = w1;
                    w2.insert(w2.begin() + pos, x);
                    g = score(d, w2, base);
                    if (keep(g)) cands.push_back({std::move(w2), g});
                }
                for (int x : shortlist(y, next, banned_y, kRobustPair)) {
                    vector<int> w2 = w1;
                    w2.insert(w2.begin() + pos + 1, x);
                    g = score(d, w2, base);
                    if (keep(g)) cands.push_back({std::move(w2), g});
                }
            }
        }
        fut.apply(d, cur_w, +1.0);
        st.add_flow(d, +1.0);
        return cands;
    }

    // Descent step on d: try the best few improving candidates.
    bool improve_demand(int d, int width, bool hot)
    {
        vector<RCand> cands = candidates(d, width, hot, [&](const GObj& g) { return worth(g); });
        if (cands.empty()) return false;
        stable_sort(cands.begin(), cands.end(), [](const RCand& a, const RCand& b) {
            if (!g_same_peak(a.g, b.g)) return g_peak_better(a.g, b.g);
            return a.g.excess < b.g.excess - kEps;
        });
        const int n = min(st.detour_try(), static_cast<int>(cands.size()));
        for (int k = 0; k < n; ++k)
            if (commit(d, cands[k].w))
                return true;
        return false;
    }

    // Install w for d without an acceptance test.
    void force(int d, const vector<int>& w)
    {
        const vector<int> old_w = st.waypoints[d];
        if (!st.apply_path(d, w)) return;
        st.recompute_obj();
        fut.apply(d, old_w, -1.0);
        fut.apply(d, w, +1.0);
        refresh();
    }

    // VNS kick: move `strength` random hot demands (warm if none) onto a
    // random candidate whose all-slot peak is at most kShakeSlack worse.
    void kick(int strength)
    {
        const vector<int>& pool = !last_hot.empty() ? last_hot : last_warm;
        if (pool.empty()) return;
        for (int k = 0; k < strength; ++k) {
            const int d = pool[st.rng() % pool.size()];
            vector<RCand> cands = candidates(d, st.shortlist_size(), true, [&](const GObj& g) {
                return g.unserved <= cur.unserved
                    && g.overflow <= cur.overflow + kShakeSlack
                    && g.peak <= cur.peak * (1.0 + kShakeSlack) + kEps;
            });
            if (!cands.empty())
                force(d, cands[st.rng() % cands.size()].w);
        }
    }

    // Slot-0 state plus the later slots' flows.
    struct Snap {
        SlotSearch::Snapshot slot0;
        vector<vector<double>> flow;
        vector<int> unserved;
        GObj g;
    };

    void save(Snap& s) const
    {
        st.save(s.slot0);
        s.flow = fut.flow;
        s.unserved = fut.unserved;
        s.g = cur;
    }

    void restore(const Snap& s)
    {
        st.restore(s.slot0);
        fut.flow = s.flow;
        fut.unserved = s.unserved;
        refresh();
    }

    // One call: demands crossing near-peak profile links first (rotating
    // window), then, if none of them moved, a window of the warm ones.
    bool step()
    {
        fut.trim();
        refresh();
        profile.assign(st.n_links, 0.0);
        for (int t = 0; t < fut.T; ++t)
            for (int i = 0; i < st.n_links; ++i)
                if (st.capacity[i] > 0.0)
                    profile[i] = max(profile[i], load_at(t, i));

        // Slots whose peak is in the near-peak band.
        thr = SlotSearch::crit_threshold(cur.peak);
        hot_slot.assign(fut.T, 0);
        for (int t = 0; t < fut.T; ++t)
            for (int i = 0; i < st.n_links && !hot_slot[t]; ++i)
                if (st.capacity[i] > 0.0 && load_at(t, i) + 1e-15 >= thr)
                    hot_slot[t] = 1;

        // Hot: the demand's route at a hot slot (its own split there)
        // crosses a near-peak link. Warm: its slot-0 route crosses a link
        // whose profile is above tau.
        vector<int> hot, warm;
        vector<double> vol(st.traffic.size(), 0.0);
        Split lk;
        for (int d = 0; d < static_cast<int>(st.traffic.size()); ++d) {
            double v = st.traffic[d];
            for (int t = 1; t < fut.T; ++t) v = max(v, fut.traffic[t][d]);
            vol[d] = v;
            if (v <= 0.0 || !st.feasible[d]) continue;
            bool is_hot = false;
            for (int t = 0; t < fut.T && !is_hot; ++t) {
                if (!hot_slot[t]) continue;
                if (t == 0) {
                    if (st.traffic[d] <= 0.0) continue;
                    lk = st.links[d];
                } else if (fut.traffic[t][d] <= 0.0 || !fut.route(t, d, st.waypoints[d], lk)) {
                    continue;
                }
                for (const auto& e : lk)
                    if (load_at(t, e.first) + 1e-15 >= thr) { is_hot = true; break; }
            }
            double m = 0.0;
            for (const auto& e : st.links[d]) m = max(m, profile[e.first]);
            if (is_hot) hot.push_back(d);
            else if (m >= tau) warm.push_back(d);
        }
        auto by_volume = [&](int a, int b) {
            if (fabs(vol[a] - vol[b]) > 1e-12) return vol[a] > vol[b];
            return a < b;
        };
        sort(hot.begin(), hot.end(), by_volume);
        sort(warm.begin(), warm.end(), by_volume);
        warm_size = static_cast<int>(warm.size());
        last_hot = hot;
        last_warm = warm;

        auto sweep = [&](const vector<int>& ds, int& cursor, int width, bool is_hot) {
            if (ds.empty()) return false;
            const int n_ds = static_cast<int>(ds.size());
            const int n = min(kExtendDemands, n_ds);
            const int start = cursor % n_ds;
            cursor = (start + n) % n_ds;
            bool any = false;
            for (int k = 0; k < n; ++k)
                if (improve_demand(ds[(start + k) % n_ds], width, is_hot))
                    any = true;
            return any;
        };
        if (sweep(hot, hot_cursor, kRobustWide, true)) return true;
        return sweep(warm, warm_cursor, st.shortlist_size(), false);
    }
};

// Run the robust phase on a finished slot-0 search state (router must hold
// the slot-0 topology). Descent until a full sweep of the hot and warm
// demands brings nothing, then VNS: keep the best local optimum, kick from
// it, descend again. Stops on time or after kRobustKicks kicks in a row
// without a new best.
static void robust_slot0(const Instance& inst, const Router& router,
                         SlotSearch& st, double time_limit, const char* tag)
{
    using clock = chrono::steady_clock;
    const auto t0 = clock::now();

    FutureSlots fut;
    fut.init(inst, router, st);
    RobustSearch rs(st, fut);
    rs.refresh();
    const double peak_before = rs.cur.peak;
    const double s0_before = st.obj.max_load;

    // Like SlotSearch::run_op, but the effect is the all-slot peak.
    auto timed = [&](int op, auto f) {
        st.begin_op(op);
        const auto s0 = clock::now();
        const double before = rs.cur.peak;
        const bool ok = f();
        if (st.stats) {
            OpStat& s = (*st.stats)[op_label(op)];
            ++s.calls;
            s.improved += ok ? 1 : 0;
            s.seconds += chrono::duration<double>(clock::now() - s0).count();
            s.peak_drop += before - rs.cur.peak;
        }
        return ok;
    };
    const bool kicks_on = !op_disabled("kick");

    RobustSearch::Snap best;
    rs.save(best);
    int idle = 0;
    int kicks = 0;
    long long calls = 0;
    while (chrono::duration<double>(clock::now() - t0).count() < time_limit) {
        ++st.iter_index;
        ++calls;
        if (timed(OP_ROBUST, [&] { return rs.step(); })) {
            idle = 0;
            continue;
        }
        const int sweeps = (rs.warm_size + kExtendDemands - 1) / kExtendDemands;
        if (++idle <= max(1, sweeps)) continue;

        // Local optimum: keep it if best, else go back to the best one.
        if (g_better(rs.cur, best.g)) {
            rs.save(best);
            kicks = 0;
        } else if (g_better(best.g, rs.cur)) {
            rs.restore(best);
        }
        if (!kicks_on || ++kicks > kRobustKicks) break;
        ++st.iter_index;
        timed(OP_KICK, [&] { rs.kick(2); return true; });
        idle = 0;
    }
    if (g_better(best.g, rs.cur)) rs.restore(best);

    cout << "  [" << tag << "] " << inst.name << " robust slot 0: "
         << "all-slot peak " << peak_before << " -> " << rs.cur.peak
         << "  slot-0 max_load " << s0_before << " -> " << st.obj.max_load
         << "  (" << calls << " calls, "
         << chrono::duration<double>(clock::now() - t0).count() << " s)\n";
}

// Copy the search state into the output Solution for a given slot.
static void fill_solution(Solution& sol, const Instance& inst, int slot,
                          SlotSearch& st)
{
    for (size_t d = 0; d < inst.demands.size(); ++d) {
        DemandPath dp;
        dp.demand_index = static_cast<int>(d);
        dp.time_slot = slot;
        dp.waypoints = st.waypoints[d];
        // Every link the demand's ECMP split touches.
        for (const auto& e : st.links[d])
            if (e.second > 1e-12) dp.links.push_back(e.first);
        dp.routed_demand = st.traffic[d];
        dp.feasible = st.feasible[d] != 0;
        dp.total_metric = 0.0;
        sol.paths.push_back(std::move(dp));
    }

    // Record per-link loads for this slot and update global max.
    for (const auto& l : inst.net.links) {
        LinkLoad ll;
        ll.time_slot = slot;
        ll.link_id = l.id;
        ll.capacity = l.capacity;
        ll.total_demand = (l.id < static_cast<int>(st.flow.size())) ? st.flow[l.id] : 0.0;
        ll.load = (l.capacity > 0.0) ? ll.total_demand / l.capacity : INF;
        if (ll.load > sol.max_link_load) {
            sol.max_link_load = ll.load;
            sol.max_link_id = l.id;
            sol.max_link_slot = slot;
        }
        sol.link_loads.push_back(ll);
    }
}

// Per-slot search budget, by instance size.
struct SlotBudget {
    int max_iters;
    double time_limit;   // seconds
};

static SlotBudget slot_budget(const Instance& inst)
{
    const size_t D = inst.demands.size();
    const size_t N = inst.net.nodes.size();
    if (D > 2000 || N > 200) return {8000, 6.0};
    if (D > 400 || N > 80) return {10000, 4.5};
    return {12000, 3.0};
}

// Shared multi-slot driver. For each slot: shortest distances under that
// slot's failures, warm start from the previous slot's waypoints, then
// `search(st, budget)` (the solver-specific part). Slot 0 is followed by the
// robust phase. Records paths, loads and the per-slot histories.
template <class SlotFn>
static Solution solve_by_slots(const Instance& inst, const char* tag, SlotFn search)
{
    using clock = chrono::steady_clock;
    const auto t_start = clock::now();

    Solution sol;
    sol.instance_name = inst.name;
    const int T = inst.num_time_slots;
    sol.max_load_history.resize(T);
    sol.max_load_history_op.resize(T);
    sol.nbhd_switches.resize(T);

    Router router;
    router.build(inst.net);
    const SlotBudget budget = slot_budget(inst);
    vector<vector<int>> prev_w;   // previous slot's waypoints

    for (int slot = 0; slot < T; ++slot) {
        unordered_set<int> forbidden;
        for (const auto& iv : inst.scenario.interventions)
            if (iv.t == slot)
                forbidden.insert(iv.links.begin(), iv.links.end());
        router.compute(forbidden);

        // Iteration counter restarts at 0 each slot (one plot per slot).
        MaxLoadHistory history;
        SlotSearch st;
        const unsigned seed = 2026u + static_cast<unsigned>(inst.index) * 1009u
                              + static_cast<unsigned>(slot) * 17u;
        st.init(inst, router, slot, slot > 0 ? &prev_w : nullptr, seed, &history);
        st.stats = slot == 0 ? &sol.op_stats_slot0 : &sol.op_stats_later;
        history.record(0, st.obj.max_load, OP_START);

        search(st, budget);

        // Slot 0: make the routes hold up in the later slots too.
        if (slot == 0 && T > 1 && !op_disabled("robust"))
            robust_slot0(inst, router, st, budget.time_limit, tag);

        fill_solution(sol, inst, slot, st);
        prev_w = st.waypoints;

        cout << "  [" << tag << "] " << inst.name << " slot " << slot
             << "  max_load=" << st.obj.max_load
             << "  unserved=" << st.obj.unserved
             << "  overflow=" << st.obj.overflow
             << "  budget=" << st.budget_used << "/"
             << (isfinite(st.budget_limit)
                     ? to_string(static_cast<int>(st.budget_limit))
                     : string("inf"))
             << "\n";

        sol.max_load_history[slot] = std::move(history.entries);
        sol.max_load_history_op[slot] = std::move(history.entry_ops);
        sol.nbhd_switches[slot] = std::move(history.switches);
    }

    sol.solve_time_seconds = chrono::duration<double>(clock::now() - t_start).count();
    return sol;
}

} // anonymous namespace

#endif // SR_CORE_H
