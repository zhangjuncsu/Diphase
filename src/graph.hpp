#pragma once

#include <queue>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_set>
#include <unordered_map>

using UID = std::size_t;

struct Queue {
    UID nid;
    int d;
    Queue() {}
    Queue(UID nid, int d): nid(nid), d(d) {}
};

// template <typename T>
class UniqueQueue {
public:
    UniqueQueue() {}
    
    bool Empty() const { return items_.empty(); }
    std::size_t Size() const { return items_.size(); }

    Queue Pop() {
        if(!Empty()) {
            Queue item = items_.front();
            items_.pop();
            ids_.erase(Val(item));
            return item;
        } else {
            std::cerr << "UniqueQueue is empty\n";
            exit(EXIT_FAILURE);
        }
    }

    void Push(const Queue &item) {
        auto val = Val(item);
        if(ids_.find(val) == ids_.end()) {
            items_.push(item);
            ids_.emplace(val);
        }
    }

    void Clear() {
        std::queue<Queue>().swap(items_);
        ids_.clear();
    }

private:
    std::size_t Val(const Queue &item) const { return item.nid << 1 | item.d; }

    std::queue<Queue> items_;
    std::unordered_set<std::size_t> ids_;
};

struct Option {
    std::string gfname;
    std::string ofname;

    int threads { 1 };
    int Check();
};

class NameMapping {
public:
    NameMapping() {}

    UID Insert(const std::string &name);
    UID TryInsert(const std::string &name);
    UID DoInsert(const std::string &name);

    std::string GetName(UID id) const;
    UID GetId(const std::string &name) const;

    bool exists(const std::string &name) const;
    bool exists(UID id) const;

    std::size_t Size() const { return names_.size(); }

private:
    std::vector<std::string> names_;
    std::unordered_map<std::string, UID> ids_;
};

class Node {
public:
    Node() {}
    Node(UID nid, std::string &seq): nid_(nid), seq_(seq) {}

    std::vector<std::size_t> GetLeftIn() const { return left_in_; }
    std::vector<std::size_t> GetRightIn() const { return right_in_; }
    std::vector<std::size_t> GetLeftOut() const { return left_out_; }
    std::vector<std::size_t> GetRightOut() const { return right_out_; }

    void AddLeftIn(std::size_t eid) { left_in_.push_back(eid); }
    void AddRightIn(std::size_t eid) { right_in_.push_back(eid); }
    void AddLeftOut(std::size_t eid) { left_out_.push_back(eid); }
    void AddRightOut(std::size_t eid) { right_out_.push_back(eid); }

    UID GetNid() const { return nid_; }
    std::string GetSeq() const { return seq_; }
    Node* GetReverse() const { return reverse_; }
    std::string GetOptional() const { return optional_; }

    bool HasOptional() const { return !optional_.empty(); }
    bool IsVisited() const { return visited_; }
    void SetVisited(bool visited) { visited_ = visited; }

    int GetPhase() const { return phase_; }
    void SetPhase(int phase) { phase_ = phase; }
    void SetReverse(Node* node) { reverse_ = node; }
    void SetSeq(std::string &seq) { seq_ = seq; }
    void SetNid(UID nid) { nid_ = nid; }
    void SetOptional(const std::string &optional) { optional_ = optional; }

    bool operator==(const Node *node) const { return nid_ == node->nid_; }

private:
    UID nid_;
    int phase_ { rand() & 1 };
    std::string seq_;
    std::string optional_;
    bool visited_ { false };
    Node* reverse_ { nullptr };
    std::vector<std::size_t> left_in_;
    std::vector<std::size_t> right_in_;
    std::vector<std::size_t> left_out_;
    std::vector<std::size_t> right_out_;
};

class ContactNode {
public:
    ContactNode() {}
    ContactNode(UID nid, int phase) {}

    UID GetNid() const { return nid_; }
    int GetPhase() const { return phase_; }

    void SetNid(UID nid) { nid_ = nid; }
    void SetPhase(int phase) { phase_ = phase; }

    std::vector<UID> GetAlts() const { return alts_; }
    void AddAlt(UID alt) { alts_.push_back(alt); }

    void Flip() { phase_ = 1 - phase_; }

private:
    UID nid_;
    int phase_;
    std::vector<UID> alts_;
};

class Edge {
public:
    Edge(UID src, char ss, UID dst, char sd, std::string &cigar): 
        src_(src), strand_src_(ss), dst_(dst), strand_dst_(sd), cigar_(cigar) {}

    Edge(UID src, char ss, UID dst, char sd, std::string &cigar, std::string &optional): 
        src_(src), strand_src_(ss), dst_(dst), strand_dst_(sd), cigar_(cigar), optional_(optional) {}

    Edge(UID src, char ss, UID dst, char sd): 
        src_(src), strand_src_(ss), dst_(dst), strand_dst_(sd) {}

    UID GetSrc() const { return src_; }
    UID GetDst() const { return dst_; }
    char GetStrandSrc() const { return strand_src_; }
    char GetStrandDst() const { return strand_dst_; }

    bool HasCigar() const { return !cigar_.empty(); }
    std::string GetCigar() const { return cigar_; }

    bool HasOptional() const { return !optional_.empty(); }
    std::string GetOptional() const { return optional_; }
    void SetOptional(const std::string &optional) { optional_ = optional; }

private:
    char strand_src_;
    char strand_dst_;
    UID src_;
    UID dst_;
    std::string cigar_;
    std::string optional_;
};

class Bubble {
public:
    Bubble() {}
    Bubble(UID from, UID to, std::vector<UID> &inside): from_(from), to_(to), inside_(inside) {}

    UID GetFrom() const { return from_; }
    UID GetTo() const { return to_; }
    std::vector<UID> GetInside() const { return inside_; }

    bool Empty() const { return inside_.empty(); }

    bool operator==(const Bubble &bubble) const { return from_ == bubble.from_ && to_ == bubble.to_; }

    friend std::ostream &operator<<(std::ostream &out, const Bubble &bubble) {
        out << "from: " << bubble.from_ << " to: " << bubble.to_ << " inside: ";
        for(auto &id: bubble.inside_) {
            out << id << " ";
        }
        out << "\n";
        return out;
    }

private:
    UID from_;
    UID to_;
    std::vector<UID> inside_;
};

class BubbleChain {
public:
    BubbleChain() {}

    std::size_t Size() const { return bubbles_.size(); }

    void AddBubble(Bubble &bubble) { bubbles_.push_back(bubble); }
    std::pair<UID, UID> FindEnds();

private:
    UID start_;
    UID end_;
    std::vector<Bubble> bubbles_;
};

class Graph {
public:
    Graph() {}
    Graph(Option &opt): opt_(opt) {}

    void AddNode(Node &node);
    void AddNode(std::vector<std::string> &items, NameMapping &nm);
    void AddEdge(std::vector<std::string> &items, NameMapping &nm);
    void ConstructGraphFromGfa(std::string &gfilename, NameMapping &nm);
    void WriteGfa(std::string &gfilename, NameMapping &nm);

    int GetNodeInDirection(UID iid, UID oid);

    Bubble FindBubble(Node &node, int d);
    std::vector<Bubble> FindBubbles();
    std::vector<BubbleChain> FindBubbleChains(std::vector<Bubble> &bubbles);

    bool HasNode(UID nid) const { return nodes_.find(nid) != nodes_.end(); }

private:
    Option opt_;
    std::unordered_map<UID, Node> nodes_;
    // std::vector<Node*> nodes_;
    std::vector<Edge> edges_;
    // std::vector<std::string> names_;
    // std::unordered_map<std::string, UID> node_name_to_id_;
};