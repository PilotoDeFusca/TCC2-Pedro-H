#include "loader.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

using json = nlohmann::json;
using namespace std;

// ---------- JSON conversion functions (must be in the GLOBAL namespace
//            so ADL can find them for our types) ----------

void from_json(const json& j, Node& n) {
    j.at("id").get_to(n.id);
    j.at("name").get_to(n.name);
}

void from_json(const json& j, Link& l) {
    j.at("id").get_to(l.id);
    j.at("from").get_to(l.from);
    j.at("to").get_to(l.to);
    j.at("metric").get_to(l.metric);
    j.at("capacity").get_to(l.capacity);
}

void from_json(const json& j, Demand& d) {
    j.at("v").get_to(d.v);
    j.at("s").get_to(d.s);
    j.at("t").get_to(d.t);
}

void from_json(const json& j, Budget& b) {
    j.at("t").get_to(b.t);
    j.at("value").get_to(b.value);
}

void from_json(const json& j, Intervention& i) {
    j.at("t").get_to(i.t);
    j.at("links").get_to(i.links);
}

void from_json(const json& j, Scenario& s) {
    j.at("max_segments").get_to(s.max_segments);
    j.at("budget").get_to(s.budget);
    j.at("interventions").get_to(s.interventions);
}

void from_json(const json& j, Network& n) {
    j.at("directed").get_to(n.directed);
    j.at("multigraph").get_to(n.multigraph);
    j.at("nodes").get_to(n.nodes);
    j.at("links").get_to(n.links);
}

// ---------- Helper: load one instance from a base path ----------
static Instance load_single_instance(const string& base_path,
                                     const string& set_name,
                                     int index)
{
    Instance inst;
    inst.set_name = set_name;
    inst.index    = index;
    inst.name     = base_path;

    // Network file
    {
        ifstream in(base_path + "-net.json");
        if (!in.is_open())
            throw runtime_error("Cannot open " + base_path + "-net.json");
        json j; in >> j;
        inst.net = j.get<Network>();
    }

    // Scenario file
    {
        ifstream in(base_path + "-scenario.json");
        if (!in.is_open())
            throw runtime_error("Cannot open " + base_path + "-scenario.json");
        json j; in >> j;
        inst.scenario = j.get<Scenario>();
    }

    // Traffic matrix file
    {
        ifstream in(base_path + "-tm.json");
        if (!in.is_open())
            throw runtime_error("Cannot open " + base_path + "-tm.json");
        json j; in >> j;
        inst.num_time_slots = j.at("num_time_slots").get<int>();
        inst.demands        = j.at("demands").get<vector<Demand>>();
    }

    return inst;
}

// ---------- Main loading function ----------
vector<Instance> load_all_instances(const vector<string>& sets,
                                    int num_instances)
{
    vector<Instance> all_instances;
    all_instances.reserve(sets.size() * num_instances);

    for (const string& set_name : sets) {
        for (int i = 1; i <= num_instances; ++i) {
            stringstream ss;
            ss << set_name << "/" << set_name << "-"
               << setw(2) << setfill('0') << i;
            string base_path = ss.str();

            try {
                Instance inst = load_single_instance(base_path, set_name, i);
                all_instances.push_back(std::move(inst));
                cout << "[OK]  Loaded " << base_path << "\n";
            } catch (const exception& e) {
                cerr << "[ERR] " << e.what() << "\n";
            }
        }
    }

    cout << "\nTotal instances loaded: " << all_instances.size() << "\n";
    return all_instances;
}