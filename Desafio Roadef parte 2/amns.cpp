#include "loader.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <unordered_set>
#include <utility>

using json = nlohmann::json;
using namespace std;

namespace {

constexpr double INF = numeric_limits<double>::infinity();

struct Arc {
    int to;
    int link_id;
    double metric;
};

struct Router {
    int num_nodes = 0;
    int max_link_id = -1;
    vector<int> node_ids;
    vector<vector<Arc>> adj;
    vector<int> pred_node; // src * num_nodes + v
    vector<int> pred_link;
    vector<double> dist;

    int idx(int src, int v) const { return src * num_nodes + v; }

    void build(const Network& net)
    {
        int max_id = -1;
        for (const auto& n : net.nodes)
            if (n.id > max_id) max_id = n.id;
        num_nodes = max_id + 1;

        node_ids.clear();
        node_ids.reserve(net.nodes.size());
        for (const auto& n : net.nodes)
            node_ids.push_back(n.id);

        adj.assign(num_nodes, {});
        max_link_id = -1;
        for (const auto& l : net.links) {
            if (l.id > max_link_id) max_link_id = l.id;
            adj[l.from].push_back({l.to, l.id, l.metric});
            if (!net.directed)
                adj[l.to].push_back({l.from, l.id, l.metric});
        }
    }

    void compute(const unordered_set<int>& forbidden)
    {
        const int N = num_nodes;
        pred_node.assign(N * N, -1);
        pred_link.assign(N * N, -1);
        dist.assign(N * N, INF);

        for (int s : node_ids) {
            vector<char> visited(N, 0);
            using P = pair<double, int>;
            priority_queue<P, vector<P>, greater<P>> pq;
            dist[idx(s, s)] = 0.0;
            pq.push({0.0, s});

            while (!pq.empty()) {
                auto [d, u] = pq.top();
                pq.pop();
                if (visited[u]) continue;
                visited[u] = 1;
                for (const auto& arc : adj[u]) {
                    if (forbidden.count(arc.link_id)) continue;
                    const double nd = d + arc.metric;
                    const int j = idx(s, arc.to);
                    if (nd < dist[j]) {
                        dist[j] = nd;
                        pred_node[j] = u;
                        pred_link[j] = arc.link_id;
                        pq.push({nd, arc.to});
                    }
                }
            }
        }
    }

    bool expand(int s, int t, const vector<int>& waypoints,
                vector<int>& links, double& metric) const
    {
        links.clear();
        metric = 0.0;
        if (s == t && waypoints.empty())
            return true;

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
            if (a < 0 || b < 0 || a >= num_nodes || b >= num_nodes)
                return false;
            if (dist[idx(a, b)] == INF)
                return false;

            vector<int> piece;
            int cur = b;
            while (cur != a) {
                const int j = idx(a, cur);
                const int link = pred_link[j];
                if (link < 0) return false;
                piece.push_back(link);
                cur = pred_node[j];
            }
            reverse(piece.begin(), piece.end());
            links.insert(links.end(), piece.begin(), piece.end());
            metric += dist[idx(a, b)];
        }
        return true;
    }
};

static uint64_t pack_seg(int u, int v)
{
    return (static_cast<uint64_t>(static_cast<uint32_t>(u)) << 32)
           | static_cast<uint32_t>(v);
}

static void segment_pairs(int s, int t, const vector<int>& w,
                          unordered_set<uint64_t>& segs)
{
    segs.clear();
    if (s == t && w.empty())
        return;
    int prev = s;
    for (int x : w) {
        if (x == prev) continue;
        segs.insert(pack_seg(prev, x));
        prev = x;
    }
    if (prev != t)
        segs.insert(pack_seg(prev, t));
}

static int segment_distance(int s, int t,
                            const vector<int>& a,
                            const vector<int>& b)
{
    unordered_set<uint64_t> sa, sb;
    segment_pairs(s, t, a, sa);
    segment_pairs(s, t, b, sb);
    int dist = 0;
    for (auto x : sa)
        if (!sb.count(x)) ++dist;
    for (auto x : sb)
        if (!sa.count(x)) ++dist;
    return dist;
}

static double budget_for_slot(const Instance& inst, int slot)
{
    if (slot <= 0)
        return numeric_limits<double>::infinity();
    for (const auto& b : inst.scenario.budget)
        if (b.t == slot)
            return b.value;
    return 0.0;
}

struct Obj {
    int unserved = 0;
    double overflow = 0.0;
    double max_load = 0.0;
    double sum_sq = 0.0;
};

static bool better(const Obj& a, const Obj& b)
{
    if (a.unserved != b.unserved) return a.unserved < b.unserved;
    if (a.overflow != b.overflow) return a.overflow < b.overflow;
    if (fabs(a.max_load - b.max_load) > 1e-12)
        return a.max_load < b.max_load;
    return a.sum_sq < b.sum_sq - 1e-12;
}

struct SlotSearch {
    const Instance* inst = nullptr;
    const Router* router = nullptr;
    int slot = 0;
    int max_waypoints = 0;
    int n_links = 0;
    vector<double> capacity;
    vector<double> flow;
    vector<double> traffic;
    vector<vector<int>> waypoints;
    vector<vector<int>> prev_waypoints;
    vector<vector<int>> links;
    vector<char> feasible;
    vector<int> budget_cost; // vs previous slot
    int budget_used = 0;
    double budget_limit = INF;
    Obj obj;
    mt19937 rng;

    void init(const Instance& in, const Router& rt, int sl,
              const vector<vector<int>>* prev, unsigned seed)
    {
        inst = &in;
        router = &rt;
        slot = sl;
        rng.seed(seed);
        max_waypoints = max(0, in.scenario.max_segments - 1);
        n_links = rt.max_link_id + 1;
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

        for (int d = 0; d < D; ++d) {
            const Demand& dem = in.demands[d];
            traffic[d] = (sl < static_cast<int>(dem.v.size())) ? dem.v[sl] : 0.0;
            if (prev)
                waypoints[d] = (*prev)[d];
            apply_path(d, waypoints[d], /*commit_flow=*/true, /*count_budget=*/true);
        }
        recompute_obj();
    }

    bool apply_path(int d, const vector<int>& w, bool commit_flow, bool count_budget)
    {
        const Demand& dem = inst->demands[d];
        vector<int> new_links;
        double metric = 0.0;
        const bool ok = router->expand(dem.s, dem.t, w, new_links, metric);

        if (count_budget) {
            const int nc = segment_distance(dem.s, dem.t, prev_waypoints[d], w);
            const int delta = nc - budget_cost[d];
            if (budget_used + delta > budget_limit + 1e-9)
                return false;
            budget_used += delta;
            budget_cost[d] = nc;
        }

        if (commit_flow && feasible[d]) {
            for (int lid : links[d])
                flow[lid] -= traffic[d];
        }

        waypoints[d] = w;
        links[d] = std::move(new_links);
        feasible[d] = ok ? 1 : 0;

        if (commit_flow && ok) {
            for (int lid : links[d])
                flow[lid] += traffic[d];
        }
        return true;
    }

    void recompute_obj()
    {
        obj.unserved = 0;
        obj.overflow = 0.0;
        obj.max_load = 0.0;
        obj.sum_sq = 0.0;
        for (size_t d = 0; d < feasible.size(); ++d) {
            if (!feasible[d] && traffic[d] > 0.0)
                ++obj.unserved;
        }
        for (int i = 0; i < n_links; ++i) {
            if (capacity[i] <= 0.0) continue;
            const double load = flow[i] / capacity[i];
            obj.max_load = max(obj.max_load, load);
            if (load > 1.0)
                obj.overflow += load - 1.0;
            obj.sum_sq += load * load;
        }
    }

    int hottest_link() const
    {
        int best = -1;
        double best_load = -1.0;
        for (int i = 0; i < n_links; ++i) {
            if (capacity[i] <= 0.0) continue;
            const double load = flow[i] / capacity[i];
            if (load > best_load) {
                best_load = load;
                best = i;
            }
        }
        return best;
    }

    int pick_demand(const vector<int>& pool)
    {
        if (pool.empty()) return -1;
        uniform_int_distribution<int> dist(0, static_cast<int>(pool.size()) - 1);
        int hot = hottest_link();
        vector<int> focused;
        if (hot >= 0) {
            for (int d : pool) {
                if (!feasible[d]) {
                    focused.push_back(d);
                    continue;
                }
                for (int lid : links[d]) {
                    if (lid == hot) {
                        focused.push_back(d);
                        break;
                    }
                }
            }
        }
        if (!focused.empty() && (rng() % 2 == 0)) {
            uniform_int_distribution<int> fd(0, static_cast<int>(focused.size()) - 1);
            return focused[fd(rng)];
        }
        return pool[dist(rng)];
    }

    int random_node_except(const unordered_set<int>& banned)
    {
        vector<int> cand;
        cand.reserve(router->node_ids.size());
        for (int id : router->node_ids)
            if (!banned.count(id))
                cand.push_back(id);
        if (cand.empty()) return -1;
        uniform_int_distribution<int> dist(0, static_cast<int>(cand.size()) - 1);
        return cand[dist(rng)];
    }

    bool neigh_change()
    {
        vector<int> pool;
        for (int d = 0; d < static_cast<int>(waypoints.size()); ++d)
            if (!waypoints[d].empty())
                pool.push_back(d);
        const int d = pick_demand(pool);
        if (d < 0) return false;

        const Demand& dem = inst->demands[d];
        uniform_int_distribution<int> pos(0, static_cast<int>(waypoints[d].size()) - 1);
        const int i = pos(rng);

        unordered_set<int> banned;
        banned.insert(dem.s);
        banned.insert(dem.t);
        banned.insert(waypoints[d][i]);
        if (i > 0) banned.insert(waypoints[d][i - 1]);
        if (i + 1 < static_cast<int>(waypoints[d].size()))
            banned.insert(waypoints[d][i + 1]);

        const int nw = random_node_except(banned);
        if (nw < 0) return false;

        vector<int> cand = waypoints[d];
        cand[i] = nw;
        return try_move(d, cand);
    }

    bool neigh_insert()
    {
        vector<int> pool;
        for (int d = 0; d < static_cast<int>(waypoints.size()); ++d)
            if (static_cast<int>(waypoints[d].size()) < max_waypoints)
                pool.push_back(d);
        const int d = pick_demand(pool);
        if (d < 0) return false;

        const Demand& dem = inst->demands[d];
        unordered_set<int> banned;
        banned.insert(dem.s);
        banned.insert(dem.t);
        for (int w : waypoints[d]) banned.insert(w);

        const int nw = random_node_except(banned);
        if (nw < 0) return false;

        vector<int> cand = waypoints[d];
        uniform_int_distribution<int> pos(0, static_cast<int>(cand.size()));
        cand.insert(cand.begin() + pos(rng), nw);
        return try_move(d, cand);
    }

    bool neigh_delete()
    {
        vector<int> pool;
        for (int d = 0; d < static_cast<int>(waypoints.size()); ++d)
            if (!waypoints[d].empty())
                pool.push_back(d);
        const int d = pick_demand(pool);
        if (d < 0) return false;

        vector<int> cand = waypoints[d];
        uniform_int_distribution<int> pos(0, static_cast<int>(cand.size()) - 1);
        cand.erase(cand.begin() + pos(rng));
        return try_move(d, cand);
    }

    bool try_move(int d, const vector<int>& cand)
    {
        const vector<int> old_w = waypoints[d];
        const vector<int> old_links = links[d];
        const char old_feas = feasible[d];
        const int old_bcost = budget_cost[d];
        const int old_bused = budget_used;
        const Obj old_obj = obj;

        if (!apply_path(d, cand, true, true)) {
            return false;
        }

        // If routing failed, still allowed if it does not worsen unserved vs previous.
        recompute_obj();
        if (better(obj, old_obj))
            return true;

        // Revert
        if (feasible[d]) {
            for (int lid : links[d])
                flow[lid] -= traffic[d];
        }
        waypoints[d] = old_w;
        links[d] = old_links;
        feasible[d] = old_feas;
        budget_cost[d] = old_bcost;
        budget_used = old_bused;
        if (old_feas) {
            for (int lid : links[d])
                flow[lid] += traffic[d];
        }
        obj = old_obj;
        return false;
    }
};

static void fill_solution(Solution& sol, const Instance& inst, int slot,
                          SlotSearch& st)
{
    for (size_t d = 0; d < inst.demands.size(); ++d) {
        DemandPath dp;
        dp.demand_index = static_cast<int>(d);
        dp.time_slot = slot;
        dp.waypoints = st.waypoints[d];
        dp.links = st.links[d];
        dp.routed_demand = st.traffic[d];
        dp.feasible = st.feasible[d] != 0;
        dp.total_metric = 0.0;
        if (dp.feasible) {
            for (int lid : dp.links) {
                // metric is already encoded in router distances; leave 0 here
                (void)lid;
            }
        }
        sol.paths.push_back(std::move(dp));
    }

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

} // anonymous namespace

Solution solve_amns(const Instance& inst)
{
    Solution sol;
    sol.instance_name = inst.name;
    sol.max_link_load = 0.0;
    sol.max_link_id = -1;
    sol.max_link_slot = -1;
    sol.total_metric = 0.0;

    Router router;
    router.build(inst.net);

    vector<vector<int>> prev_w;
    const bool has_prev = false;
    (void)has_prev;

    using clock = chrono::steady_clock;

    for (int slot = 0; slot < inst.num_time_slots; ++slot) {
        unordered_set<int> forbidden;
        for (const auto& iv : inst.scenario.interventions) {
            if (iv.t == slot)
                forbidden.insert(iv.links.begin(), iv.links.end());
        }
        router.compute(forbidden);

        SlotSearch st;
        const unsigned seed = 2026u + static_cast<unsigned>(inst.index) * 1009u
                              + static_cast<unsigned>(slot) * 17u;
        st.init(inst, router, slot, slot > 0 ? &prev_w : nullptr, seed);

        const int D = static_cast<int>(inst.demands.size());
        const int N = static_cast<int>(inst.net.nodes.size());
        int max_iters = 12000;
        double time_limit = 3.0;
        if (D > 2000 || N > 200) {
            max_iters = 8000;
            time_limit = 6.0;
        } else if (D > 400 || N > 80) {
            max_iters = 10000;
            time_limit = 4.5;
        }

        array<double, 3> score = {1.0, 1.0, 1.0};
        const auto t0 = clock::now();

        for (int it = 0; it < max_iters; ++it) {
            const double elapsed = chrono::duration<double>(clock::now() - t0).count();
            if (elapsed > time_limit) break;

            const double sum = score[0] + score[1] + score[2];
            uniform_real_distribution<double> ur(0.0, sum);
            double r = ur(st.rng);
            int k = 0;
            if (r < score[0]) k = 0;
            else if (r < score[0] + score[1]) k = 1;
            else k = 2;

            bool improved = false;
            for (int attempt = 0; attempt < 3 && !improved; ++attempt) {
                if (k == 0) improved = st.neigh_change();
                else if (k == 1) improved = st.neigh_insert();
                else improved = st.neigh_delete();
                if (!improved)
                    k = (k + 1) % 3;
            }

            for (double& s : score) s *= 0.995;
            if (improved)
                score[k] += 1.0;
            else
                score[k] = max(0.2, score[k] * 0.97);
        }

        fill_solution(sol, inst, slot, st);
        prev_w = st.waypoints;

        cout << "  [AMNS] " << inst.name << " slot " << slot
             << "  max_load=" << st.obj.max_load
             << "  unserved=" << st.obj.unserved
             << "  overflow=" << st.obj.overflow
             << "  budget=" << st.budget_used << "/"
             << (isfinite(st.budget_limit) ? to_string(static_cast<int>(st.budget_limit)) : string("inf"))
             << "\n";
    }

    return sol;
}

void write_srpaths(const Solution& sol, const string& path)
{
    json j;
    j["srpaths"] = json::array();

    for (const auto& dp : sol.paths) {
        json e;
        e["d"] = dp.demand_index;
        e["t"] = dp.time_slot;
        e["w"] = dp.waypoints;
        j["srpaths"].push_back(std::move(e));
    }

    ofstream out(path);
    if (!out.is_open()) {
        cerr << "Cannot write " << path << "\n";
        return;
    }
    out << j.dump(2) << "\n";
}
