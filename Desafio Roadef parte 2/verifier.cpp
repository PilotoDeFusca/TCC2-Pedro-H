#include "loader.h"

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
    int prev = s;
    for (int x : w) {
        if (x == prev) continue;
        segs.insert(pack_seg(prev, x));
        prev = x;
    }
    if (prev != t)
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

// For each slot, get the path assigned to demand d (or nullptr if none)

VerificationReport verify_solution(const Solution& sol, const Instance& inst)
{
    VerificationReport rep;
    rep.instance_name = inst.name;
    rep.ok = true;

    // --- 1. Capacity violations ---
    // For each (slot, link), check load > 1.0
    for (const auto& ll : sol.link_loads) {
        if (ll.load > 1.0 + 1e-9) {
            rep.ok = false;
            VerificationIssue iss;
            iss.kind         = "capacity";
            iss.time_slot    = ll.time_slot;
            iss.demand_index = -1;
            iss.link_id      = ll.link_id;
            ostringstream oss;
            oss << "Link " << ll.link_id
                << " at slot " << ll.time_slot
                << " has load " << ll.load
                << " (" << ll.total_demand << "/" << ll.capacity << ")";
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

    if (rep.ok) {
        cout << "  OK: no issues detected.\n";
        return;
    }

    cout << "  FAILED: " << rep.issues.size() << " issue(s)\n";

    // Group by kind for readability
    for (const string& kind : {"capacity", "unserved", "intervention", "segments", "budget"}) {
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