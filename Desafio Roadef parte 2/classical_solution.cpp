// classical_solution.cpp
// Classical Variable Neighborhood Descent (VND) for multi-slot Segment
// Routing path optimization.
//
// Difference vs. solve_amns():
//   - AMNS  : samples one move type per iteration, adaptively weighted.
//   - VND   : insert, then change, then safe-delete (restart at insert after
//             any improvement), each restricted to the critical cut with
//             delta-scored shortlists. At a local optimum, VNS shakes from the
//             best solution seen and descends again.
//
// Everything else (ECMP Router, SlotSearch state, budget, objective,
// apply_path, try_move, robust slot 0, fill_solution) is shared with amns.cpp
// through sr_core.h.

#include "sr_core.h"

using namespace std;

namespace {

// ---------------------------------------------------------------------------
// VND / VNS driver for a single slot
// ---------------------------------------------------------------------------
//
// Neighborhoods in descent order (insert first: on slot 0 every route is the
// IGP path, so change/delete have nothing to work on until insert has run):
//   insert (main attack on the bottleneck)
//   insert2 (two waypoints at once, when one cannot leave the cut)
//   extend (slot 0 only: lengthen routes through cold links)
//   change (diversify an existing detour)
//   delete (safe cleanup / budget recovery)
//
// Classic VND: after any improvement, restart from k = 0; stop at a local
// optimum of all three. Then VNS: keep the best local optimum, shake from it
// with increasing strength, and descend again until kMaxStall shakes in a
// row bring no improvement (or time/iterations run out).

static void vnd_descent_slot(SlotSearch& st,
                             int max_iters_total, double time_limit_sec)
{
    using clock = chrono::steady_clock;
    const auto t0 = clock::now();
    auto out_of_budget = [&](int calls) {
        return calls >= max_iters_total ||
               chrono::duration<double>(clock::now() - t0).count() > time_limit_sec;
    };

    constexpr int kMaxStall = 8;
    SlotSearch::Snapshot best;
    st.save(best);
    int calls = 0;
    int stall = 0;
    bool stop = false;

    // Extend has no segment budget to spend only at slot 0. Moves switched
    // off with SR_DISABLE_OPS are left out.
    vector<int> ops;
    for (int op : (st.slot == 0)
             ? vector<int>{OP_INSERT, OP_INSERT2, OP_EXTEND, OP_CHANGE, OP_DELETE}
             : vector<int>{OP_INSERT, OP_INSERT2, OP_CHANGE, OP_DELETE})
        if (!op_disabled(op_label(op))) ops.push_back(op);
    const int n_ops = static_cast<int>(ops.size());
    const bool shake_on = !op_disabled("shake");

    while (!stop) {
        // ---- VND descent ----
        int k = 0;
        while (k < n_ops) {
            if (out_of_budget(calls)) { stop = true; break; }
            ++st.iter_index;
            ++calls;
            const int op = ops[k];
            const bool improved = st.run_op(op, [&] { return st.neighborhood(op); });
            k = improved ? 0 : k + 1;
        }

        // ---- move or not ----
        if (better(st.obj, best.obj)) {
            st.save(best);
            stall = 0;
        } else {
            if (better(best.obj, st.obj)) st.restore(best);
            ++stall;
        }
        if (stop || !shake_on || stall >= kMaxStall || out_of_budget(calls)) break;

        // ---- shake from the best basin ----
        ++st.iter_index;
        ++calls;
        st.run_op(OP_SHAKE, [&] { st.shake(2 + stall); return true; });
    }

    if (better(best.obj, st.obj)) st.restore(best);
}

} // anonymous namespace

// Entry point: run the VND across all time slots and return the solution.
Solution solve_classical(const Instance& inst)
{
    return solve_by_slots(inst, "CLASSICAL", [](SlotSearch& st, const SlotBudget& b) {
        vnd_descent_slot(st, b.max_iters, b.time_limit);
    });
}
