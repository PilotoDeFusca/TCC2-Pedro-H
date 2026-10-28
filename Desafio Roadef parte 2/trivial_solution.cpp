// trivial_solution.cpp
// Baseline: a network without segment-routing optimization.
//
// Every demand has an empty waypoint list in every slot, so its traffic
// follows plain IGP routing: ECMP over the shortest paths of that slot's
// topology (failed links removed, so traffic re-routes around them). No
// segment budget is used. Links can end up above capacity (load > 1); those
// are kept as they are and reported as viability problems by the verifier.
//
// Loads are computed with the same ECMP routing as the solvers (sr_core.h),
// so they match the checker.

#include "sr_core.h"

using namespace std;

Solution solve_shortest_paths(const Instance& inst)
{
    using clock = chrono::steady_clock;
    const auto t_start = clock::now();
    const int T = inst.num_time_slots;

    Router router;
    router.build(inst.net);
    const vector<vector<int>> no_waypoints(inst.demands.size());

    Solution sol;
    sol.instance_name = inst.name;
    for (int slot = 0; slot < T; ++slot) {
        unordered_set<int> forbidden;
        for (const auto& iv : inst.scenario.interventions)
            if (iv.t == slot)
                forbidden.insert(iv.links.begin(), iv.links.end());
        router.compute(forbidden);

        SlotSearch st;
        st.init(inst, router, slot, &no_waypoints, 0, nullptr);
        fill_solution(sol, inst, slot, st);

        int over = 0;
        for (const auto& ll : sol.link_loads)
            if (ll.time_slot == slot && ll.load > 1.0 + 1e-9) ++over;
        cout << "  [TRIVIAL] " << inst.name << " slot " << slot
             << "  max_load=" << st.obj.max_load
             << "  unserved=" << st.obj.unserved
             << "  links_over_capacity=" << over << "\n";
    }
    sol.max_load_history.assign(T, {});
    sol.max_load_history_op.assign(T, {});
    sol.nbhd_switches.assign(T, {});
    sol.solve_time_seconds = chrono::duration<double>(clock::now() - t_start).count();
    return sol;
}

void print_solution(const Solution& sol, const Instance& inst)
{
    int feasible_count = 0;
    int infeasible_count = 0;

    cout << "--- Demand paths (time slot 0) ---\n";
    for (const auto& dp : sol.paths) {
        const Demand& dem = inst.demands[dp.demand_index];

        cout << "Demand #" << dp.demand_index
             << " (s=" << dem.s << ", t=" << dem.t << ")"
             << "  demand=" << dp.routed_demand
             << "  -> ";

        if (!dp.feasible) {
            cout << "NO PATH FOUND\n";
            ++infeasible_count;
            continue;
        }

        ++feasible_count;

        cout << "metric=" << dp.total_metric
             << "  path=[";

        for (size_t k = 0; k < dp.links.size(); ++k) {
            if (k) cout << " -> ";
            cout << dp.links[k];
        }
        cout << "]\n";
    }

    cout << "\n--- Link loads ---\n";
    cout << "link_id  capacity   demand      load\n";
    for (const auto& ll : sol.link_loads) {
        cout << "  " << ll.link_id
             << "      " << ll.capacity
             << "      " << ll.total_demand
             << "      " << ll.load << "\n";
    }

    cout << "\nSummary for " << sol.instance_name << ": "
         << feasible_count << " feasible, "
         << infeasible_count << " infeasible"
         << " (out of " << sol.paths.size() << " demand-time pairs)\n";

    cout << "\n--- Objective ---\n";
    cout << "Max link load = " << sol.max_link_load
         << "  (link " << sol.max_link_id
         << ", slot " << sol.max_link_slot << ")\n";
}
