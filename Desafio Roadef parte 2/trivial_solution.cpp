#include "loader.h"

#include <iostream>
#include <queue>
#include <limits>
#include <algorithm>
#include <unordered_map>

using namespace std;

namespace {

constexpr double INF = std::numeric_limits<double>::infinity();

// A simple adjacency-list entry: (neighbor_node_id, link_id, metric)
struct Arc {
    int to;
    int link_id;
    double metric;
};

// Build adjacency list from the network (directed or undirected)
vector<vector<Arc>> build_adjacency(const Network& net)
{
    int max_id = -1;
    for (const auto& n : net.nodes)
        if (n.id > max_id) max_id = n.id;

    vector<vector<Arc>> adj(max_id + 1);

    for (const auto& l : net.links) {
        adj[l.from].push_back({l.to, l.id, l.metric});
        if (!net.directed) {
            adj[l.to].push_back({l.from, l.id, l.metric});
        }
    }
    return adj;
}

// Reconstruct the link-id path from predecessor maps
vector<int> reconstruct_path(int s, int t,
                             const vector<int>& pred_node,
                             const vector<int>& pred_link)
{
    vector<int> path;
    int cur = t;
    while (cur != s) {
        int link = pred_link[cur];
        if (link < 0) return {};
        path.push_back(link);
        cur = pred_node[cur];
    }
    reverse(path.begin(), path.end());
    return path;
}

// Dijkstra from source s to target t. Returns the sequence of link ids.
vector<int> dijkstra(const vector<vector<Arc>>& adj,
                     int num_nodes,
                     int s, int t,
                     double& out_metric)
{
    if (s < 0 || t < 0 ||
        s >= static_cast<int>(adj.size()) ||
        t >= static_cast<int>(adj.size())) {
        out_metric = INF;
        return {};
    }

    vector<double> dist(num_nodes, INF);
    vector<int>    pred_node(num_nodes, -1);
    vector<int>    pred_link(num_nodes, -1);
    vector<bool>   visited(num_nodes, false);

    using P = pair<double, int>;
    priority_queue<P, vector<P>, greater<P>> pq;

    dist[s] = 0.0;
    pq.push({0.0, s});

    while (!pq.empty()) {
        P top = pq.top();
        pq.pop();
        double d = top.first;
        int u = top.second;

        if (visited[u]) continue;
        visited[u] = true;
        if (u == t) break;

        for (const auto& arc : adj[u]) {
            double nd = d + arc.metric;
            if (nd < dist[arc.to]) {
                dist[arc.to]      = nd;
                pred_node[arc.to] = u;
                pred_link[arc.to] = arc.link_id;
                pq.push(P(nd, arc.to));
            }
        }
    }

    out_metric = dist[t];
    if (dist[t] == INF) return {};
    return reconstruct_path(s, t, pred_node, pred_link);
}

} // anonymous namespace

Solution solve_shortest_paths(const Instance& inst)
{
    Solution sol;
    sol.instance_name = inst.name;
    sol.max_link_load = 0.0;
    sol.max_link_id   = -1;
    sol.max_link_slot = -1;
    sol.total_metric  = 0.0;

    int max_id = -1;
    for (const auto& n : inst.net.nodes)
        if (n.id > max_id) max_id = n.id;
    const int num_nodes = max_id + 1;

    auto adj = build_adjacency(inst.net);

    // Iterate every time slot
    for (int slot = 0; slot < inst.num_time_slots; ++slot) {

        // Fresh link-demand map per slot
        unordered_map<int, double> link_demand;
        for (const auto& l : inst.net.links)
            link_demand[l.id] = 0.0;

        // Route each demand for this slot
        for (size_t d = 0; d < inst.demands.size(); ++d) {
            const Demand& dem = inst.demands[d];

            DemandPath dp;
            dp.demand_index = static_cast<int>(d);
            dp.time_slot    = slot;
            dp.waypoints    = {};

            double traffic = (slot < static_cast<int>(dem.v.size()))
                             ? dem.v[slot]
                             : 0.0;
            dp.routed_demand = traffic;

            double metric = INF;
            dp.links = dijkstra(adj, num_nodes, dem.s, dem.t, metric);
            dp.total_metric = metric;
            dp.feasible     = !dp.links.empty() || (dem.s == dem.t);

            if (dp.feasible) {
                sol.total_metric += metric;
                for (int lid : dp.links)
                    link_demand[lid] += traffic;
            }

            sol.paths.push_back(std::move(dp));
        }

        // Build LinkLoad entries for this slot
        for (const auto& l : inst.net.links) {
            LinkLoad ll;
            ll.time_slot    = slot;
            ll.link_id      = l.id;
            ll.capacity     = l.capacity;
            ll.total_demand = link_demand[l.id];
            ll.load = (l.capacity > 0.0)
                      ? ll.total_demand / l.capacity
                      : INF;

            if (ll.load > sol.max_link_load) {
                sol.max_link_load = ll.load;
                sol.max_link_id   = l.id;
                sol.max_link_slot = slot;
            }

            sol.link_loads.push_back(ll);
        }
    }

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

    cout << "\n--- Objective ---\n";
    cout << "Max link load = " << sol.max_link_load
         << "  (link " << sol.max_link_id << ")\n";
    cout << "Total metric  = " << sol.total_metric << "\n";

    cout << "\nSummary for " << sol.instance_name << ": "
         << feasible_count << " feasible, "
         << infeasible_count << " infeasible"
         << " (out of " << sol.paths.size() << " demand-time pairs)\n";

    cout << "\n--- Objective ---\n";
    cout << "Max link load = " << sol.max_link_load
          << "  (link " << sol.max_link_id
          << ", slot " << sol.max_link_slot << ")\n";
    cout << "Total metric  = " << sol.total_metric << "\n";
}