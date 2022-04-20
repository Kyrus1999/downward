#ifndef HEURISTICS_ADDITIVE_HEURISTIC_H
#define HEURISTICS_ADDITIVE_HEURISTIC_H

#include "relaxation_heuristic.h"

#include "../algorithms/priority_queues.h"
#include "../utils/collections.h"

#include <cassert>

class State;

namespace additive_heuristic {
using relaxation_heuristic::PropID;
using relaxation_heuristic::OpID;

using relaxation_heuristic::NO_OP;

using relaxation_heuristic::Proposition;
using relaxation_heuristic::UnaryOperator;
using relaxation_heuristic::EffectNode;
using relaxation_heuristic::InnerNode;

class AdditiveHeuristic : public relaxation_heuristic::RelaxationHeuristic {
    /* Costs larger than MAX_COST_VALUE are clamped to max_value. The
       precise value (100M) is a bit of a hack, since other parts of
       the code don't reliably check against overflow as of this
       writing. With a value of 100M, we want to ensure that even
       weighted A* with a weight of 10 will have f values comfortably
       below the signed 32-bit int upper bound.
     */
    static const int MAX_COST_VALUE = 100000000;

    priority_queues::AdaptiveQueue<PropID> queue;
    bool did_write_overflow_warning;
    void (*enque_function)(relaxation_heuristic::RelaxationHeuristic*, PropID, int, OpID);

    void setup_exploration_queue(const State &state);
    void relaxed_exploration();
    void mark_preferred_operators(const State &state, PropID goal_id);

    void enqueue_if_necessary(PropID prop_id, int cost, OpID op_id) {
        if (prop_id == relaxation_heuristic::NO_PROP)
            return;
        assert(cost >= 0);
        Proposition *prop = get_proposition(prop_id);
        if (prop->cost == -1 || prop->cost > cost) {
            prop->cost = cost;
            prop->reached_by = op_id;
            queue.push(cost, prop_id);
        }
        assert(prop->cost != -1 && prop->cost <= cost);
    }

    static void enque(relaxation_heuristic::RelaxationHeuristic *additive_heuristic, PropID prop_id, int cost, OpID op_id){
        ((AdditiveHeuristic*) additive_heuristic)->enqueue_if_necessary(prop_id, cost, op_id);
    }

    void increase_cost(int &cost, int amount) {
        assert(cost >= 0);
        assert(amount >= 0);
        cost += amount;
        if (cost > MAX_COST_VALUE) {
            write_overflow_warning();
            cost = MAX_COST_VALUE;
        }
    }

    void write_overflow_warning();
protected:
    virtual int compute_heuristic(const State &ancestor_state) override;

    // Common part of h^add and h^ff computation.
    int compute_add_and_ff(const State &state);
public:
    explicit AdditiveHeuristic(const options::Options &opts);

    /*
      TODO: The two methods below are temporarily needed for the CEGAR
      heuristic. In the long run it might be better to split the
      computation from the heuristic class. Then the CEGAR code could
      use the computation object instead of the heuristic.
    */
    void compute_heuristic_for_cegar(const State &state);

    int get_cost_for_cegar(int var, int value) const {
        return get_proposition(var, value)->cost;
    }
};
}

#endif
