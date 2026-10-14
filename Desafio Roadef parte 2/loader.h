#ifndef LOADER_H
#define LOADER_H

#include <vector>
#include <string>

#include <json.hpp>

// ---------- Data structures ----------
struct Node {
    int id;
    std::string name;
};

struct Link {
    int id;
    int from;
    int to;
    double metric;
    double capacity;
};

struct Demand {
    std::vector<double> v;
    int s;
    int t;
};

struct Budget {
    int t;
    double value;
};

struct Intervention {
    int t;
    std::vector<int> links;
};

struct Scenario {
    int max_segments;
    std::vector<Budget> budget;
    std::vector<Intervention> interventions;
};

struct Network {
    bool directed;
    bool multigraph;
    std::vector<Node> nodes;
    std::vector<Link> links;
};

struct Instance {
    std::string name;
    std::string set_name;
    int         index;
    Network     net;
    Scenario    scenario;
    int         num_time_slots;
    std::vector<Demand> demands;
};

// ---------- Solution representation ----------
struct DemandPath {
    int demand_index;
    int time_slot;
    std::vector<int> waypoints;  // SR nodes between source and destination
    std::vector<int> links;      // expanded IGP link sequence
    double total_metric;
    double routed_demand;
    bool   feasible;
};

struct LinkLoad {
    int    time_slot;      // NEW
    int    link_id;
    double capacity;
    double load;
    double total_demand;
};

struct Solution {
    std::string instance_name;
    double max_link_load = 0.0;
    int    max_link_id   = -1;
    int    max_link_slot = -1;
    double total_metric  = 0.0;
    double solve_time_seconds = 0.0;

    // Per-slot max-load history.
    // max_load_history[slot] = list of (iteration_within_slot, max_load)
    // recorded whenever the accepted solution's max link load dropped.
    std::vector<std::vector<std::pair<long long, double>>> max_load_history;

    std::vector<DemandPath> paths;
    std::vector<LinkLoad>   link_loads;
};

// ---------- Loading function ----------
std::vector<Instance> load_all_instances(
    const std::vector<std::string>& sets = {"setA", "setB"},
    int num_instances = 20);

// ---------- Algorithm ----------
Solution solve_shortest_paths(const Instance& inst);
Solution solve_amns(const Instance& inst);
Solution solve_classical(const Instance& inst);
void     print_solution(const Solution& sol, const Instance& inst);
void     write_srpaths(const Solution& sol, const std::string& path);
void     write_solution(const Solution& sol, const std::string& path);

// ---------- Verifier ----------
struct VerificationIssue {
    std::string kind;        // "capacity", "unserved", "intervention", "segments"
    int time_slot;
    int demand_index;        // -1 if not applicable
    int link_id;             // -1 if not applicable
    std::string message;
};

struct VerificationReport {
    std::string instance_name;
    bool ok;                         // true = no issues at all
    std::vector<VerificationIssue> issues;
};

VerificationReport verify_solution(const Solution& sol, const Instance& inst);
void               print_verification(const VerificationReport& rep);

#endif // LOADER_H