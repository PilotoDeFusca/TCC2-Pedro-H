#include "loader.h"

#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

static void print_slot_loads(const Solution& sol, const Instance& inst)
{
    vector<double> slot_max(inst.num_time_slots, 0.0);
    vector<int> slot_link(inst.num_time_slots, -1);

    for (const auto& ll : sol.link_loads) {
        if (ll.time_slot < 0 || ll.time_slot >= inst.num_time_slots)
            continue;
        if (ll.load > slot_max[ll.time_slot]) {
            slot_max[ll.time_slot] = ll.load;
            slot_link[ll.time_slot] = ll.link_id;
        }
    }

    cout << "\n--- Maximum load per time slot ---\n";
    cout << fixed << setprecision(6);
    for (int t = 0; t < inst.num_time_slots; ++t) {
        cout << "  " << inst.name << "  slot " << t
             << "  max_load=" << slot_max[t]
             << "  (link " << slot_link[t] << ")\n";
    }
}

static void print_sr_nodes(const Solution& sol, const Instance& inst)
{
    cout << "\n--- Segment routing nodes (non-empty) ---\n";
    int shown = 0;
    const int cap = 40;
    for (const auto& dp : sol.paths) {
        if (dp.waypoints.empty()) continue;
        if (shown >= cap) {
            cout << "  ... (further SR lists written to " << inst.name
                 << "-srpaths.json)\n";
            break;
        }
        cout << "  slot=" << dp.time_slot
             << " demand=" << dp.demand_index
             << " w=[";
        for (size_t i = 0; i < dp.waypoints.size(); ++i) {
            if (i) cout << ", ";
            cout << dp.waypoints[i];
        }
        cout << "]\n";
        ++shown;
    }
    if (shown == 0)
        cout << "  (all demands use empty waypoint lists: IGP shortest path)\n";
}

int main() {
    vector<Instance> instances = load_all_instances({"setA","setB"}, 20);

    cout << "\n==================== FINAL SUMMARY ====================\n";
    vector<VerificationReport> reports;

    for (const auto& inst : instances) {
        cout << "\n=====================================================\n";
        cout << "Solving " << inst.name
             << "  (nodes=" << inst.net.nodes.size()
             << ", links=" << inst.net.links.size()
             << ", demands=" << inst.demands.size()
             << ", slots=" << inst.num_time_slots
             << ", max_segments=" << inst.scenario.max_segments << ")\n";
        cout << "=====================================================\n";

        Solution sol = solve_amns(inst);
        write_srpaths(sol, inst.name + "-srpaths.json");
        print_sr_nodes(sol, inst);
        print_slot_loads(sol, inst);

        cout << "Instance max load = " << sol.max_link_load
             << "  (link " << sol.max_link_id
             << ", slot " << sol.max_link_slot << ")\n";

        VerificationReport rep = verify_solution(sol, inst);
        print_verification(rep);
        reports.push_back(rep);
    }

    cout << "\n=====================================================\n";
    cout << "Maximum load of each instance and time slot\n";
    cout << "=====================================================\n";

    // Re-print a compact table. Loads were already printed while solving;
    // reload is not needed because we printed per instance. Collect from files
    // is unnecessary: reprint verification status here.
    int failed = 0;
    for (const auto& rep : reports) {
        if (!rep.ok) {
            ++failed;
            cout << "  " << rep.instance_name << "  VIABLE? NO  ("
                 << rep.issues.size() << " issue(s))\n";
            for (const auto& iss : rep.issues)
                cout << "      [" << iss.kind << "] " << iss.message << "\n";
        } else {
            cout << "  " << rep.instance_name << "  VIABLE? YES\n";
        }
    }

    if (failed == 0)
        cout << "\nNo viability errors on any instance.\n";
    else
        cout << "\nViability errors on " << failed << " instance(s).\n";

    return 0;
}
