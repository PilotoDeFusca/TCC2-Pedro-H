#include "loader.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>

using namespace std;

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
    // Like the official checker, every consecutive pair is a segment, even a
    // repeated waypoint (x, x): the solvers never emit one, and this catches
    // it if they do.
    int prev = s;
    for (int x : w) {
        segs.insert(pack_seg(prev, x));
        prev = x;
    }
    segs.insert(pack_seg(prev, t));
}

static int segment_distance(int s, int t, const vector<int>& a, const vector<int>& b)
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

// -------- Helpers --------

// For each slot, get the path assigned to demand d (or nullptr if none)
static const DemandPath* find_path(const Solution& sol, int slot, int demand_index)
{
    for (const auto& dp : sol.paths) {
        if (dp.time_slot == slot && dp.demand_index == demand_index)
            return &dp;
    }
    return nullptr;
}

// Traffic per link that no routing can avoid at `slot`: from a source whose
// only usable out-link is that link (or into a target whose only usable
// in-link is), followed along single-link chains. Same rule as the solvers.
static unordered_map<int, double> forced_flow(const Instance& inst, int slot)
{
    unordered_set<int> down;
    for (const auto& iv : inst.scenario.interventions)
        if (iv.t == slot) down.insert(iv.links.begin(), iv.links.end());

    int N = 0;
    for (const auto& n : inst.net.nodes) N = max(N, n.id + 1);
    vector<vector<pair<int, int>>> out(N), in(N);   // (link id, other end)
    for (const auto& l : inst.net.links) {
        if (down.count(l.id)) continue;
        out[l.from].push_back({l.id, l.to});
        in[l.to].push_back({l.id, l.from});
        if (!inst.net.directed) {
            out[l.to].push_back({l.id, l.from});
            in[l.from].push_back({l.id, l.to});
        }
    }

    unordered_map<int, double> forced;
    for (const auto& dem : inst.demands) {
        const double traffic = (slot < (int)dem.v.size()) ? dem.v[slot] : 0.0;
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

VerificationReport verify_solution(const Solution& sol, const Instance& inst)
{
    VerificationReport rep;
    rep.instance_name = inst.name;
    rep.ok = true;

    // --- 1. Capacity violations ---
    // For each (slot, link), check load > 1.0. If every unit of flow on the
    // link is forced (no routing could avoid it), report it as
    // "capacity_unavoidable": a note that does not make the solution fail.
    unordered_map<int, unordered_map<int, double>> forced_by_slot;
    for (const auto& ll : sol.link_loads) {
        if (ll.load > 1.0 + 1e-9) {
            auto fit = forced_by_slot.find(ll.time_slot);
            if (fit == forced_by_slot.end())
                fit = forced_by_slot.emplace(ll.time_slot,
                                             forced_flow(inst, ll.time_slot)).first;
            auto lit = fit->second.find(ll.link_id);
            const double forced = (lit != fit->second.end()) ? lit->second : 0.0;
            const bool unavoidable = ll.total_demand <= forced * (1.0 + 1e-9) + 1e-12;

            if (!unavoidable)
                rep.ok = false;
            VerificationIssue iss;
            iss.kind         = unavoidable ? "capacity_unavoidable" : "capacity";
            iss.time_slot    = ll.time_slot;
            iss.demand_index = -1;
            iss.link_id      = ll.link_id;
            ostringstream oss;
            oss << "Link " << ll.link_id
                << " at slot " << ll.time_slot
                << " has load " << setprecision(10) << ll.load << setprecision(6)
                << " (" << ll.total_demand << "/" << ll.capacity << ")";
            if (unavoidable)
                oss << ", all of it forced (only usable link for its demands)";
            else if (forced > 0.0)
                oss << ", forced lower bound " << forced / ll.capacity;
            iss.message = oss.str();
            rep.issues.push_back(iss);
        }
    }

    // --- 2. Unserved demands ---
    // For each demand and each slot, check there is a feasible path.
    // A zero-traffic demand counts as served (nothing to route).
    for (int slot = 0; slot < inst.num_time_slots; ++slot) {
        for (size_t d = 0; d < inst.demands.size(); ++d) {
            const Demand& dem = inst.demands[d];
            double traffic = (slot < (int)dem.v.size()) ? dem.v[slot] : 0.0;

            const DemandPath* dp = find_path(sol, slot, (int)d);

            if (dp == nullptr) {
                rep.ok = false;
                VerificationIssue iss;
                iss.kind         = "unserved";
                iss.time_slot    = slot;
                iss.demand_index = (int)d;
                iss.link_id      = -1;
                iss.message = "No path entry for demand";
                rep.issues.push_back(iss);
                continue;
            }

            if (!dp->feasible && traffic > 0.0) {
                rep.ok = false;
                VerificationIssue iss;
                iss.kind         = "unserved";
                iss.time_slot    = slot;
                iss.demand_index = (int)d;
                iss.link_id      = -1;
                ostringstream oss;
                oss << "Demand #" << d
                    << " (s=" << dem.s << ", t=" << dem.t << ")"
                    << " at slot " << slot
                    << " has traffic " << traffic
                    << " but no feasible path";
                iss.message = oss.str();
                rep.issues.push_back(iss);
            }
        }
    }

    // --- 3. Forbidden links (intervention links) ---
    // Build the set of forbidden link ids *per slot* from the scenario.
    // The scenario lists interventions with a time index "t". We interpret
    // t as the time slot during which the listed links are *intervened*
    // (i.e. unavailable for routing).
    unordered_map<int, vector<int>> forbidden; // slot -> link ids
    for (const auto& iv : inst.scenario.interventions) {
        forbidden[iv.t].insert(forbidden[iv.t].end(),
                               iv.links.begin(), iv.links.end());
    }

    for (const auto& dp : sol.paths) {
        auto it = forbidden.find(dp.time_slot);
        if (it == forbidden.end()) continue;

        const auto& fset = it->second;
        for (int lid : dp.links) {
            if (std::find(fset.begin(), fset.end(), lid) != fset.end()) {
                rep.ok = false;
                VerificationIssue iss;
                iss.kind         = "intervention";
                iss.time_slot    = dp.time_slot;
                iss.demand_index = dp.demand_index;
                iss.link_id      = lid;
                ostringstream oss;
                oss << "Demand #" << dp.demand_index
                    << " at slot " << dp.time_slot
                    << " uses forbidden link " << lid;
                iss.message = oss.str();
                rep.issues.push_back(iss);
            }
        }
    }

    // --- 4. SR segments vs max_segments ---
    // Number of segments = |waypoints| + 1 (source -> first WP -> ... -> dest)
    const int max_seg = inst.scenario.max_segments;
    for (const auto& dp : sol.paths) {
        if (!dp.feasible) continue;
        const int nseg = static_cast<int>(dp.waypoints.size()) + 1;
        if (nseg > max_seg) {
            rep.ok = false;
            VerificationIssue iss;
            iss.kind         = "segments";
            iss.time_slot    = dp.time_slot;
            iss.demand_index = dp.demand_index;
            iss.link_id      = -1;
            ostringstream oss;
            oss << "Demand #" << dp.demand_index
                << " at slot " << dp.time_slot
                << " has " << nseg
                << " SR segments, exceeds max_segments=" << max_seg;
            iss.message = oss.str();
            rep.issues.push_back(iss);
        }
    }

    // --- 5. Budget between consecutive time slots ---
    unordered_map<int, double> budget_of;
    for (const auto& b : inst.scenario.budget)
        budget_of[b.t] = b.value;

    for (int slot = 1; slot < inst.num_time_slots; ++slot) {
        double limit = 0.0;
        auto bit = budget_of.find(slot);
        if (bit != budget_of.end())
            limit = bit->second;

        int cost = 0;
        for (size_t d = 0; d < inst.demands.size(); ++d) {
            const DemandPath* prev = find_path(sol, slot - 1, static_cast<int>(d));
            const DemandPath* cur  = find_path(sol, slot, static_cast<int>(d));
            if (!prev || !cur) continue;
            const Demand& dem = inst.demands[d];
            cost += segment_distance(dem.s, dem.t, prev->waypoints, cur->waypoints);
        }

        if (cost > limit + 1e-9) {
            rep.ok = false;
            VerificationIssue iss;
            iss.kind         = "budget";
            iss.time_slot    = slot;
            iss.demand_index = -1;
            iss.link_id      = -1;
            ostringstream oss;
            oss << "Budget violated at slot " << slot
                << ": cost=" << cost << " > budget=" << limit;
            iss.message = oss.str();
            rep.issues.push_back(iss);
        }
    }

    return rep;
}

// -------- Pretty print --------

void print_verification(const VerificationReport& rep)
{
    cout << "\n=== Verification of " << rep.instance_name << " ===\n";

    if (rep.ok && rep.issues.empty()) {
        cout << "  OK: no issues detected.\n";
        return;
    }

    if (rep.ok)
        cout << "  OK: " << rep.issues.size() << " note(s), no avoidable issues\n";
    else
        cout << "  FAILED: " << rep.issues.size() << " issue(s)\n";

    // Group by kind for readability
    for (const string& kind : {"capacity", "capacity_unavoidable", "unserved",
                               "intervention", "segments", "budget"}) {
        int n = 0;
        for (const auto& iss : rep.issues)
            if (iss.kind == kind) ++n;
        if (n == 0) continue;

        cout << "\n  [" << kind << "] " << n << " issue(s):\n";
        for (const auto& iss : rep.issues) {
            if (iss.kind != kind) continue;
            cout << "    - " << iss.message << "\n";
        }
    }
}