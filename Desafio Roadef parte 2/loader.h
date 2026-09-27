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
    std::vector<int> links;
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
    std::vector<DemandPath> paths;

    std::vector<LinkLoad> link_loads;  // now one entry per (slot, link)
    double max_link_load;              // worst over all slots
    int    max_link_id;
    int    max_link_slot;              // NEW
    double total_metric;
};

// ---------- Loading function ----------
std::vector<Instance> load_all_instances(
    const std::vector<std::string>& sets = {"setA", "setB"},
    int num_instances = 20);

// ---------- Algorithm ----------
Solution solve_shortest_paths(const Instance& inst);
void     print_solution(const Solution& sol, const Instance& inst);

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