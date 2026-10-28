// amns.cpp
// Adaptive Multi-Neighborhood Search (AMNS) heuristic for multi-slot
// Segment-Routing (SR) path optimization.
//
// Goal: for each time slot, choose a small set of waypoints per demand so
// that link congestion is minimized, while respecting a per-slot "budget"
// on how many segments may change relative to the previous slot.
//
// Routing follows the checker: between consecutive waypoints, traffic is
// split by ECMP over all shortest paths (see Router in sr_core.h).
//
// High-level pipeline per slot:
//   1. Compute all-pairs shortest distances under current failed links.
//   2. Warm-start waypoints from the previous slot.
//   3. Local search restricted to the critical cut (links with
//      load >= max_load - eps, plus unserved demands), largest volume first.
//        insert : main attack on the bottleneck
//        change : diversification when insert stalls
//        delete : safe cleanup / segment-budget recovery
//      Each demand is taken off the network, a residual-ranked shortlist of
//      detour nodes is delta-scored exactly (new max load, avoids hot links,
//      budget cost), and only the best few are tried with try_move.
//   4. Accept strict improvements of (unserved, overflow, max_load,
//      #critical links, sum_sq). Plateau moves on the peak are allowed via
//      #critical links / sum_sq; later slots do not spend budget on sum_sq
//      alone. Equal-objective moves that free budget are accepted.
//   5. On stall: VNS shake (slightly worsening moves allowed) from the best
//      solution seen; the best solution is restored at the end of the slot.
//   6. Slot 0 only: robust phase (robust_slot0 in sr_core.h). Budgets are too
//      small to re-route much later, so slot-0 routes are re-scored against
//      every slot's traffic and failures and moved onto low-load nodes.
//   7. Record paths + link loads, carry waypoints to next slot.

#include "sr_core.h"

#include <filesystem>
#include <fstream>
#include <thread>

using json = nlohmann::json;
using namespace std;

namespace {

// AMNS stops a slot once this many neighborhood calls pass without progress
// on the peak (see peak_progress). Later slots, with their small budgets,
// usually stall within a few dozen calls.
constexpr int kPatience = 100;

// Progress for the stopping rule: fewer unserved demands, less overflow, a
// lower max load, or the same peak on fewer links. Tie-breaker-only gains
// (warm-link excess, sum of squares) do not keep the search alive.
static bool peak_progress(const Obj& n, const Obj& ref)
{
    if (peak_better(n, ref)) return true;
    return same_peak(n, ref) && n.n_crit < ref.n_crit;
}

// AMNS on one slot: each iteration samples a neighborhood (insert, change,
// delete, extend, insert2) by adaptive weight, falling back to the others if
// it fails; after two failures in a row, VNS shake from the best solution
// seen.
static void amns_slot(SlotSearch& st, const SlotBudget& budget)
{
    using clock = chrono::steady_clock;
    const auto t0 = clock::now();

    // Neighborhoods and their initial adaptive weights. Extend only runs at
    // slot 0, where there is no segment budget. Moves switched off with
    // SR_DISABLE_OPS are left out.
    vector<int> op_id;
    vector<double> score;
    const vector<pair<int, double>> all_ops = {
        {OP_INSERT, 2.0}, {OP_CHANGE, 1.0}, {OP_DELETE, 0.55},
        {OP_EXTEND, 1.0}, {OP_INSERT2, 1.0}};
    for (const auto& [op, w] : all_ops) {
        if (op == OP_EXTEND && st.slot != 0) continue;
        if (op_disabled(op_label(op))) continue;
        op_id.push_back(op);
        score.push_back(w);
    }
    const int n_ops = static_cast<int>(op_id.size());
    if (n_ops == 0) return;
    const bool shake_on = !op_disabled("shake");
    int fail_streak = 0;
    SlotSearch::Snapshot best;
    st.save(best);
    Obj progress_ref = st.obj;             // best peak reached so far
    long long last_progress = st.iter_index;

    for (int it = 0; it < budget.max_iters; ++it) {
        if (chrono::duration<double>(clock::now() - t0).count() > budget.time_limit)
            break;
        if (st.iter_index - last_progress > kPatience)
            break;

        double sum = 0.0;
        for (int o = 0; o < n_ops; ++o) sum += score[o];
        uniform_real_distribution<double> ur(0.0, sum);
        double r = ur(st.rng);
        int k = 0;
        while (k + 1 < n_ops && r >= score[k]) {
            r -= score[k];
            ++k;
        }

        bool improved = false;
        for (int attempt = 0; attempt < n_ops && !improved; ++attempt) {
            // One iteration = one neighborhood call, so fallbacks to
            // another operator show up as their own iteration.
            ++st.iter_index;
            improved = st.run_op(op_id[k], [&] { return st.neighborhood(op_id[k]); });
            if (!improved)
                k = (k + 1) % n_ops;
        }

        if (peak_progress(st.obj, progress_ref)) {
            progress_ref = st.obj;
            last_progress = st.iter_index;
        }

        for (double& s : score) s *= 0.995;
        if (improved) {
            score[k] += 1.0;
            fail_streak = 0;
        } else {
            score[k] = max(0.2, score[k] * 0.97);
            if (++fail_streak >= 2) {
                // VNS: remember the best local optimum and always
                // perturb from it, never from a worse basin.
                if (better(st.obj, best.obj)) {
                    st.save(best);
                } else if (better(best.obj, st.obj)) {
                    st.restore(best);
                }
                if (shake_on) {
                    ++st.iter_index;
                    st.run_op(OP_SHAKE, [&] { st.shake(2); return true; });
                }
                fail_streak = 0;
            }
        }
    }

    // Return the best solution seen in this slot, not the last one.
    if (better(best.obj, st.obj)) st.restore(best);
}

} // anonymous namespace

// Entry point: run AMNS across all time slots and return the solution.
Solution solve_amns(const Instance& inst)
{
    return solve_by_slots(inst, "AMNS", amns_slot);
}

// Write `text` to `path` through a temporary file and a rename, retrying a
// few times: in a OneDrive-synced folder an existing file can be briefly
// locked, and opening it directly then fails. Reports failure on stderr.
static void write_file(const string& path, const string& text)
{
    namespace fs = std::filesystem;
    const string tmp = path + ".tmp";
    for (int attempt = 0; attempt < 8; ++attempt) {
        {
            ofstream out(tmp, ios::binary | ios::trunc);
            if (out.is_open()) out << text;
        }
        std::error_code ec;
        fs::rename(tmp, path, ec);
        if (!ec) return;
        this_thread::sleep_for(chrono::milliseconds(250 * (attempt + 1)));
    }
    std::error_code ec;
    fs::remove(tmp, ec);
    cerr << "ERROR: cannot write " << path << "\n";
}

void write_solution(const Solution& sol, const string& path)
{
    json j;
    j["instance_name"]      = sol.instance_name;
    j["max_link_load"]      = sol.max_link_load;
    j["max_link_id"]        = sol.max_link_id;
    j["max_link_slot"]      = sol.max_link_slot;
    j["total_metric"]       = sol.total_metric;
    j["solve_time_seconds"] = sol.solve_time_seconds;

    // Per-demand paths.
    j["paths"] = json::array();
    for (const auto& dp : sol.paths) {
        json e;
        e["d"]             = dp.demand_index;
        e["t"]             = dp.time_slot;
        e["w"]             = dp.waypoints;
        e["links"]         = dp.links;
        e["routed_demand"] = dp.routed_demand;
        e["feasible"]      = dp.feasible;
        e["total_metric"]  = dp.total_metric;
        j["paths"].push_back(std::move(e));
    }

    // Per-(slot, link) loads.
    j["link_loads"] = json::array();
    for (const auto& ll : sol.link_loads) {
        json e;
        e["t"]            = ll.time_slot;
        e["link_id"]      = ll.link_id;
        e["capacity"]     = ll.capacity;
        e["total_demand"] = ll.total_demand;
        e["load"]         = ll.load;
        j["link_loads"].push_back(std::move(e));
    }

    // Per-slot max-load history:
    //   max_load_history[slot] = array of {iter, max_load, op}
    //   nbhd_switches[slot]    = array of {iter, op}  (op active from iter on)
    // op is one of "start", "insert", "change", "delete", "shake", "extend",
    // "robust", "kick".
    const auto op_name = op_label;
    j["max_load_history"] = json::array();
    for (size_t s = 0; s < sol.max_load_history.size(); ++s) {
        json slot_arr = json::array();
        const auto& ops = (s < sol.max_load_history_op.size())
                              ? sol.max_load_history_op[s] : std::vector<int>{};
        for (size_t i = 0; i < sol.max_load_history[s].size(); ++i) {
            const auto& [it, val] = sol.max_load_history[s][i];
            json e;
            e["iter"]     = it;
            e["max_load"] = val;
            e["op"]       = op_name(i < ops.size() ? ops[i] : OP_START);
            slot_arr.push_back(std::move(e));
        }
        j["max_load_history"].push_back(std::move(slot_arr));
    }
    j["nbhd_switches"] = json::array();
    for (const auto& slot_sw : sol.nbhd_switches) {
        json slot_arr = json::array();
        for (const auto& [it, op] : slot_sw) {
            json e;
            e["iter"] = it;
            e["op"]   = op_name(op);
            slot_arr.push_back(std::move(e));
        }
        j["nbhd_switches"].push_back(std::move(slot_arr));
    }

    // Move statistics: op_stats.{slot0,later}.<move> = {calls, improved,
    // seconds, peak_drop} (see OpStat in loader.h).
    auto stats_json = [](const map<string, OpStat>& m) {
        json o = json::object();
        for (const auto& [name, s] : m)
            o[name] = {{"calls", s.calls}, {"improved", s.improved},
                       {"seconds", s.seconds}, {"peak_drop", s.peak_drop}};
        return o;
    };
    j["op_stats"] = {{"slot0", stats_json(sol.op_stats_slot0)},
                     {"later", stats_json(sol.op_stats_later)}};

    write_file(path, j.dump(2) + "\n");
}

// Dump only the SR paths (demand index, time slot, waypoints) to JSON.
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

    write_file(path, j.dump(2) + "\n");
}