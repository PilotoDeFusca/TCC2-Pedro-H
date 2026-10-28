#include "loader.h"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

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
            cout << "  ... (further SR lists written to "
                 << inst.name << "-srpaths.json)\n";
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

// ---------------------------------------------------------------------------
// Export a solution with the requested split:
//
//   full solution JSON  -> <sol_dir>/<instance_base>-solution.json
//   SR paths JSON       -> <sr_dir>/<instance_base>-srpaths.json
//
// Example for inst.name = "setA/setA-01":
//   amns_solution/setA-01-solution.json
//   setA_amns/setA-01-srpaths.json
// ---------------------------------------------------------------------------
static void export_solution(const Solution& sol,
                            const Instance& inst,
                            const fs::path& sol_dir,
                            const fs::path& sr_dir,
                            const char* label)
{
    fs::create_directories(sol_dir);
    fs::create_directories(sr_dir);

    // Base filename: strip the "setX/" prefix from inst.name.
    string base_name = inst.name;
    const auto slash = base_name.find('/');
    if (slash != string::npos)
        base_name = base_name.substr(slash + 1);

    const fs::path sol_path  = sol_dir / (base_name + "-solution.json");
    const fs::path sr_path   = sr_dir  / (base_name + "-srpaths.json");

    write_solution(sol, sol_path.string());
    write_srpaths(sol, sr_path.string());

    cout << "  [" << label << "] wrote " << sol_path.string() << "\n";
    cout << "  [" << label << "] wrote " << sr_path.string()  << "\n";
}

static void print_report_summary(const char* label,
                                 const vector<VerificationReport>& reps)
{
    cout << "\n=====================================================\n";
    cout << "Verification summary (" << label << ")\n";
    cout << "=====================================================\n";
    int failed = 0;
    for (const auto& rep : reps) {
        if (!rep.ok) {
            ++failed;
            cout << "  " << rep.instance_name << "  VIABLE? NO  ("
                 << rep.issues.size() << " issue(s))\n";
        } else if (!rep.issues.empty()) {
            cout << "  " << rep.instance_name << "  VIABLE? YES  ("
                 << rep.issues.size() << " unavoidable overload(s))\n";
        } else {
            cout << "  " << rep.instance_name << "  VIABLE? YES\n";
        }
        for (const auto& iss : rep.issues)
            cout << "      [" << iss.kind << "] " << iss.message << "\n";
    }
    if (failed == 0)
        cout << "\nNo viability errors on any instance.\n";
    else
        cout << "\nViability errors on " << failed << " instance(s).\n";
}

int main() {
    vector<Instance> instances = load_all_instances({"setA", "setB", "setC", "setX"}, 20);

    // One entry per algorithm: full solutions go to a flat folder
    // (<sol_dir>/), SR paths to one folder per set (<set>_<sr_suffix>/).
    struct Algo {
        const char* label;          // log label
        const char* summary;        // verification summary title
        Solution (*solve)(const Instance&);
        fs::path sol_dir;
        string sr_suffix;
        vector<VerificationReport> reports;
    };
    vector<Algo> algos = {
        {"AMNS", "AMNS", solve_amns, "amns_solution", "amns", {}},
        {"VND", "Classical VND", solve_classical, "vnd_solution", "vnd", {}},
        {"TRIVIAL", "Trivial (IGP, no waypoints)", solve_shortest_paths,
         "trivial_solution", "trivial", {}},
    };

    cout << "\n==================== FINAL SUMMARY ====================\n";
    for (const auto& inst : instances) {
        cout << "\n=====================================================\n";
        cout << "Solving " << inst.name
             << "  (nodes=" << inst.net.nodes.size()
             << ", links=" << inst.net.links.size()
             << ", demands=" << inst.demands.size()
             << ", slots=" << inst.num_time_slots
             << ", max_segments=" << inst.scenario.max_segments << ")\n";
        cout << "=====================================================\n";

        for (auto& algo : algos) {
            cout << "\n>>> Running " << algo.summary << " ...\n";
            const Solution sol = algo.solve(inst);
            export_solution(sol, inst, algo.sol_dir,
                            fs::path(inst.set_name + "_" + algo.sr_suffix), algo.label);
            print_sr_nodes(sol, inst);
            print_slot_loads(sol, inst);
            cout << algo.label << " instance max load = " << sol.max_link_load
                 << "  (link " << sol.max_link_id
                 << ", slot " << sol.max_link_slot << ")\n";
            algo.reports.push_back(verify_solution(sol, inst));
        }
    }

    for (const auto& algo : algos)
        print_report_summary(algo.summary, algo.reports);
    return 0;
}

