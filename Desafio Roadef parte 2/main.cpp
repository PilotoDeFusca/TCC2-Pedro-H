#include "loader.h"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

// Forward declarations (defined in amns.cpp / classical_solution.cpp)
Solution solve_amns(const Instance& inst);
Solution solve_classical(const Instance& inst);
void write_solution(const Solution& sol, const string& path);
void write_srpaths(const Solution& sol, const string& path);

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

int main() {
    // Two sets, two instances each.  NOTE: two separate strings.
    vector<Instance> instances = load_all_instances({"setA", "setB", "setC", "setX"}, 20);

    // Two flat folders for the full solution JSONs.
    const fs::path amns_sol_dir = "amns_solution";
    const fs::path vnd_sol_dir  = "vnd_solution";

    // One srpaths folder per (set, algorithm).  Created on demand below.
    auto sr_dir_for = [](const string& set_name, const string& algo) {
        return fs::path(set_name + "_" + algo);
    };

    fs::create_directories(amns_sol_dir);
    fs::create_directories(vnd_sol_dir);

    cout << "\n==================== FINAL SUMMARY ====================\n";
    vector<VerificationReport> reports_amns;
    vector<VerificationReport> reports_vnd;

    for (const auto& inst : instances) {
        cout << "\n=====================================================\n";
        cout << "Solving " << inst.name
             << "  (nodes=" << inst.net.nodes.size()
             << ", links=" << inst.net.links.size()
             << ", demands=" << inst.demands.size()
             << ", slots=" << inst.num_time_slots
             << ", max_segments=" << inst.scenario.max_segments << ")\n";
        cout << "=====================================================\n";

        // Derive the set name ("setA", "setB") from the instance.
        // Prefer inst.set_name; fall back to the prefix of inst.name.
        string set_name = inst.set_name;
        if (set_name.empty()) {
            const auto slash = inst.name.find('/');
            set_name = (slash == string::npos)
                           ? string("unknown")
                           : inst.name.substr(0, slash);
        }

        // -----------------------------------------------------------------
        // (A) AMNS: solve ONCE, then export.
        //     Solution -> amns_solution/
        //     SR paths -> setX_amns/
        // -----------------------------------------------------------------
        cout << "\n>>> Running AMNS ...\n";
        Solution sol_amns = solve_amns(inst);
        export_solution(sol_amns,
                        inst,
                        amns_sol_dir,
                        sr_dir_for(set_name, "amns"),
                        "AMNS");

        print_sr_nodes(sol_amns, inst);
        print_slot_loads(sol_amns, inst);

        cout << "AMNS instance max load = " << sol_amns.max_link_load
             << "  (link " << sol_amns.max_link_id
             << ", slot " << sol_amns.max_link_slot << ")\n";

        reports_amns.push_back(verify_solution(sol_amns, inst));

        // -----------------------------------------------------------------
        // (B) Classical VND: solve ONCE, then export.
        //     Solution -> vnd_solution/
        //     SR paths -> setX_vnd/
        // -----------------------------------------------------------------
        cout << "\n>>> Running Classical VND ...\n";
        Solution sol_vnd = solve_classical(inst);
        export_solution(sol_vnd,
                        inst,
                        vnd_sol_dir,
                        sr_dir_for(set_name, "vnd"),
                        "VND");

        print_sr_nodes(sol_vnd, inst);
        print_slot_loads(sol_vnd, inst);

        cout << "VND instance max load = " << sol_vnd.max_link_load
             << "  (link " << sol_vnd.max_link_id
             << ", slot " << sol_vnd.max_link_slot << ")\n";

        reports_vnd.push_back(verify_solution(sol_vnd, inst));
    }

    // ---------------------------------------------------------------------
    // Verification summary for both algorithms.
    // ---------------------------------------------------------------------
    auto print_report_summary = [](const char* label,
                                   const vector<VerificationReport>& reps) {
        cout << "\n=====================================================\n";
        cout << "Verification summary (" << label << ")\n";
        cout << "=====================================================\n";
        int failed = 0;
        for (const auto& rep : reps) {
            if (!rep.ok) {
                ++failed;
                cout << "  " << rep.instance_name << "  VIABLE? NO  ("
                     << rep.issues.size() << " issue(s))\n";
                for (const auto& iss : rep.issues)
                    cout << "      [" << iss.kind << "] "
                         << iss.message << "\n";
            } else {
                cout << "  " << rep.instance_name << "  VIABLE? YES\n";
            }
        }
        if (failed == 0)
            cout << "\nNo viability errors on any instance.\n";
        else
            cout << "\nViability errors on " << failed
                 << " instance(s).\n";
    };

    print_report_summary("AMNS", reports_amns);
    print_report_summary("Classical VND", reports_vnd);

    return 0;
}