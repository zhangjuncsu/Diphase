#include "graph.hpp"
#include "file_io.hpp"
#include "utility.hpp"
#include "logger.hpp"

#include <queue>
#include <cassert>
#include <ostream>
#include <unordered_set>

#include <getopt.h>

UID ParseGfaNodeId(NameMapping &name_mapping, std::string &name) {
    if(name_mapping.exists(name)) {
        return name_mapping.GetId(name);
    } else {
        return name_mapping.Insert(name);
    }
}

UID NameMapping::Insert(const std::string &name) {
    if(exists(name)) {
        LOG(ERROR)("name %s already exists\n", name.c_str());
    }
    return DoInsert(name);
}

UID NameMapping::TryInsert(const std::string &name) {
    auto result = ids_.find(name);
    if(result == ids_.end()) {
        return DoInsert(name);
    } else {
        return result->second;
    }
}

UID NameMapping::DoInsert(const std::string &name) {
    names_.push_back(name);
    UID id = names_.size() - 1;
    ids_.emplace(name, id);
    return id;
}

std::string NameMapping::GetName(UID id) const {
    if(id >= names_.size()) {
        LOG(ERROR)("id %zd does not exist\n", id);
    }
    return names_[id];
}

UID NameMapping::GetId(const std::string &name) const {
    auto result = ids_.find(name);
    if(result == ids_.end()) {
        LOG(ERROR)("name %s does not exist\n", name.c_str());
    }
    return result->second;
}

bool NameMapping::exists(const std::string &name) const {
    return ids_.find(name) != ids_.end();
}

bool NameMapping::exists(UID id) const {
    return id < names_.size();
}

void Graph::ConstructGraphFromGfa(std::string &gfilename, NameMapping &name_mapping) {
    LOG(INFO)("constructing graph from %s\n", gfilename.c_str());
    GzFileReader reader(gfilename);
    if(reader.Valid()) {
        std::string line = reader.GetNoEmptyLine();
        while(!line.empty()) {
            auto items = SplitString(line, '\t');
            if(items[0] == "S") {
                AddNode(items, name_mapping);
            } else if(items[0] == "L") {
                AddEdge(items, name_mapping);
            }
            line = reader.GetNoEmptyLine();
        }
    } else {
        std::cerr << "[" << GetCurTime() << "] Failed to open " << gfilename << " for reading\n";
        exit(EXIT_FAILURE);
    }
}

void Graph::AddNode(std::vector<std::string> &items, NameMapping &name_mapping) {
    if(items.size() < 3) {
        LOG(ERROR)("gfa node %s contains no sequence\n", items[1].c_str());
    }
    if(items[3] == "*") {
        LOG(WARNING)("gfa node %s has no sequence\n", items[1].c_str());
        return;
    }
    UID nid = ParseGfaNodeId(name_mapping, items[1]);
    Node node(nid, items[2]);
    if(items.size() > 3) {
        std::ostringstream oss;
        for(std::size_t i = 3; i < items.size(); ++i) {
            oss << items[i] << "\t";
        }
        node.SetOptional(oss.str());
    }
    nodes_[nid] = node;
}

void Graph::AddEdge(std::vector<std::string> &items, NameMapping &name_mapping) {
    if(items.size() < 6) {
        LOG(ERROR)("error link format in gfa\n");
    }
    UID src_id = ParseGfaNodeId(name_mapping, items[1]);
    UID dst_id = ParseGfaNodeId(name_mapping, items[3]);
    auto& src = nodes_[src_id];
    auto& dst = nodes_[dst_id];
    edges_.emplace_back(src_id, items[2][0], dst_id, items[4][0], items[5]);
    if(items.size() > 6) {
        std::ostringstream oss;
        for(std::size_t i = 6; i < items.size(); ++i) {
            oss << items[i] << "\t";
        }
        edges_.back().SetOptional(oss.str());
    }
    auto edge_id = edges_.size() - 1;
    if(items[2][0] == '+') {
        src.AddRightOut(dst_id);
    } else {
        src.AddLeftOut(dst_id);
    }
    if(items[4][0] == '+') {
        dst.AddLeftIn(src_id);
    } else {
        dst.AddRightIn(src_id);
    }
}

void Graph::WriteGfa(std::string &gfilename, NameMapping &name_mapping) {
    std::ofstream out(gfilename);
    LOG(INFO)("writing graph to %s\n", gfilename.c_str());
    LOG(INFO)("number of nodes: %zd\n", nodes_.size());
    LOG(INFO)("number of edges: %zd\n", edges_.size());
    if(out.is_open()) {
        out << "H\tVN:Z:1.0\n";
        std::vector<UID> node_ids;
        for(auto node: nodes_) {
            node_ids.push_back(node.first);
        }
        std::sort(node_ids.begin(), node_ids.end());
        for(auto id: node_ids) {
            auto name = name_mapping.GetName(id);
            auto seq = nodes_[id].GetSeq();
            out << "S\t" << name << "\t" << seq;
            if(nodes_[id].HasOptional()) out << "\t" << nodes_[id].GetOptional();
            out << "\n";
            // std::cerr << "node: " << id << " " << name << " ";
            // for(auto li: nodes_[id].GetLeftIn()) {
            //     std::cerr << li << " ";
            // }
            // std::cerr << "| ";
            // for(auto ri: nodes_[id].GetRightIn()) {
            //     std::cerr << ri << " ";
            // }
            // std::cerr << "| ";
            // for(auto lo: nodes_[id].GetLeftOut()) {
            //     std::cerr << lo << " ";
            // }
            // std::cerr << "| ";
            // for(auto ro: nodes_[id].GetRightOut()) {
            //     std::cerr << ro << " ";
            // }
            // std::cerr << "\n";
        }
        for(auto edge: edges_) {
            auto src = edge.GetSrc();
            auto dst = edge.GetDst();
            auto src_name = name_mapping.GetName(src);
            auto dst_name = name_mapping.GetName(dst);
            auto ss = edge.GetStrandSrc();
            auto sd = edge.GetStrandDst();
            out << "L\t" << src_name << "\t" << ss << "\t" << dst_name << "\t" << sd << "\t" << edge.GetCigar();
            if(edge.HasOptional()) out << "\t" << edge.GetOptional();
            out << "\n";
        }
    } else {
        std::cerr << "[" << GetCurTime() << "] Failed to open " << gfilename << " for writing\n";
        exit(EXIT_FAILURE);
    }
    out.close();
}

int Graph::GetNodeInDirection(UID iid, UID oid) {
    auto left_in = nodes_[iid].GetLeftIn();
    auto right_in = nodes_[iid].GetRightIn();
    if(std::find(left_in.begin(), left_in.end(), oid) != left_in.end()) return 0;
    else if(std::find(right_in.begin(), right_in.end(), oid) != right_in.end()) return 1;
    else LOG(ERROR)("node %zd has no edge from node %zd\n", iid, oid);
}

Bubble Graph::FindBubble(Node &start, int d) {
    std::unordered_set<UID> seen;
    std::unordered_set<UID> visited;
    std::vector<UID> inside;
    seen.insert(start.GetNid());
    UniqueQueue queue;
    queue.Push(Queue(start.GetNid(), d));
    while(!queue.Empty()) {
        auto item = queue.Pop();
        visited.insert(item.nid);
        inside.push_back(item.nid);
        seen.erase(item.nid);
        std::vector<UID> children;
        if(item.d == 0) children = nodes_[item.nid].GetRightOut();
        else children = nodes_[item.nid].GetLeftOut();
        if(children.size() == 0) break;
        for(auto child: children) {
            int u_child_dir;
            std::vector<UID> u_parants;
            if(GetNodeInDirection(child, item.nid) == 0) {
                u_child_dir = 0;
                u_parants = nodes_[child].GetLeftIn();
            } else {
                u_child_dir = 1;
                u_parants = nodes_[child].GetRightIn();
            }
            if(child == start.GetNid()) {
                queue.Clear();
                break;
            }
            seen.insert(child);
            bool all_seen = true;
            for(auto u_parant: u_parants) {
                if(visited.find(u_parant) == visited.end()) {
                    all_seen = false;
                    break;
                }
            }
            if(all_seen) {
                queue.Push(Queue(child, u_child_dir));
            }
        }
        if(queue.Size() == 1 && seen.size() == 1) {
            auto item = queue.Pop();
            inside.push_back(item.nid);
            if(inside.size() == 2) break;
            inside.erase(std::find(inside.begin(), inside.end(), start.GetNid()));
            inside.erase(std::find(inside.begin(), inside.end(), item.nid));
            Bubble bubble(start.GetNid(), item.nid, inside);
            return bubble;
        }
    }
    return Bubble();
}

std::vector<Bubble> Graph::FindBubbles() {
    std::vector<Bubble> bubbles;
    for(auto& node: nodes_) {
        if(node.second.IsVisited()) continue;
        for(int d: { 0, 1 }) {
            auto bubble = FindBubble(node.second, d);
            if(!bubble.Empty()) {
                std::cerr << "bubble: " << bubble.GetFrom() << " " << bubble.GetTo() << "\n";
                for(auto id: bubble.GetInside()) {
                    std::cerr << id << " ";
                }
                std::cerr << "\n";
                bubbles.push_back(bubble);
            }
        }
    }
    return bubbles;
}

std::vector<BubbleChain> Graph::FindBubbleChains(std::vector<Bubble> &bubbles) {
    std::vector<BubbleChain> chains;
    std::unordered_map<UID, std::vector<Bubble>> bubble_map;
    for(auto bubble: bubbles) {
        bubble_map[bubble.GetFrom()].push_back(bubble);
        bubble_map[bubble.GetTo()].push_back(bubble);
    }
    std::vector<UID> start;
    for(auto& item: bubble_map) {
        if(item.second.size() == 1) {
            start.push_back(item.first);
        }
    }
    std::sort(start.begin(), start.end());
    for(auto s: start) {
        BubbleChain chain;
        Bubble current_bubble;
        if(bubble_map[s].size() != 0) current_bubble = bubble_map[s][0];
        else continue;
        UID current_node = s;

        while(true) {
            UID next_node;
            if(current_node = current_bubble.GetFrom()) next_node = current_bubble.GetTo();
            else next_node = current_bubble.GetFrom();
            bubble_map[next_node].erase(std::find(bubble_map[next_node].begin(), bubble_map[next_node].end(), current_bubble));
            chain.AddBubble(current_bubble);
            if(bubble_map[next_node].size() < 1) break;
            current_bubble = bubble_map[next_node][0];
            current_node = next_node;
        }

        if(chain.Size() != 0) {
            auto ends = chain.FindEnds();
            chains.push_back(chain);
        }
    }
    return chains;
}

int Option::Check() {
    if(gfname.empty()) {
        std::cerr << "[" << GetCurTime() << "] Graph file name is not specified\n";
        return 1;
    }
    if(threads < 1) {
        LOG(WARNING)("threads should be at least 1, set to 1\n");
        threads = 1;
    }
    return 0;
}

void USAGE(Option &opt) {
    std::cerr << "Usage: phasing [options]\n";
    std::cerr << "Options:\n";
    std::cerr << "\t\t-g | --graph   <FILE>    graph file in gfa format\n";
    std::cerr << "\t\t-t | --threads <INT>     number of threads [1]\n";
}

int ParseArgument(int argc, char **argv, Option &opt) {
    struct option long_options[] = {
        {"graph",   required_argument,  0,  'g'},
        {"threads", required_argument,  0,  't'},
        {"output",  required_argument,  0,  'o'},
        {"help",    no_argument,        0,  'h'},
        {0, 0, 0, 0}
    };
    const char* short_options = "g:t:o:h";
    int c;
    while((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1) {
        if(c == 'g') opt.gfname = optarg;
        else if(c == 't') opt.threads = atoi(optarg);
        else if(c == 'o') opt.ofname = optarg;
        else if(c == 'h') {
            USAGE(opt);
            exit(EXIT_SUCCESS);
        } else if(c == ':') {
            std::cerr << "[" << GetCurTime() << "] ERROR missing argument in option " << optopt << "\n";
            return 1;
        } else if(c == '?') {
            std::cerr << "[" << GetCurTime() << "] ERROR unknown option " << optopt << "\n";
            return 1;
        }
    }
    return opt.Check();
}

int main(int argc, char **argv) {
    Option opt;
    if(ParseArgument(argc, argv, opt)) return 1;
    NameMapping name_mapping;
    Graph graph(opt);
    graph.ConstructGraphFromGfa(opt.gfname, name_mapping);
    graph.WriteGfa(opt.ofname, name_mapping);
    auto bubbles = graph.FindBubbles();
    graph.FindBubbleChains(bubbles);
    return 0;
}