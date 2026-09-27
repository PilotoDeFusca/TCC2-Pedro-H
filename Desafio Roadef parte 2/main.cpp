#include "loader.h"

#include <iostream>

using namespace std;

int main() {
    vector<Instance> instances = load_all_instances({"setA"}, 1);

    for (const auto& inst : instances) {
        cout << "\n=====================================================\n";
        cout << "Solving " << inst.name << "\n";
        cout << "=====================================================\n";

        Solution sol = solve_shortest_paths(inst);
        print_solution(sol, inst);

        VerificationReport rep = verify_solution(sol, inst);
        print_verification(rep);
    }

    return 0;
}