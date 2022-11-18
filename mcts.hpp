/**
 * @file mcts.hpp
 * @brief Single header of the \a mcts library.
 * Implements a framework allowing to apply the Monte-Carlo Tree Search to your concrete problems.
 *
 * Basically, the MCTS class implements the "Pseudo-code for Monte-Carlo Tree Search" (\see
 * https://dke.maastrichtuniversity.nl/m.winands/documents/pMCTS.pdf) using virtual classes.
 * Therefore, as long as you provide an implementation of those classes that suits your needs,
 * everything should be fine.
 *
 * Here is the list of classes that need to be implemented :
 *<ul>
 *<li> State : Should be able to represent any state of the game.</li>
 *<li> Move : Represents a move of the game. When performed (on a given \a State), it alters the
 *State of the game.</li>
 *<li> TerminalCriteria : This allow MCTS to determinate a final state. </li>
 *<li> TerminalEval : This allow MCTS to evaluate a final state (and backpropagate the eval). <:li>
 *<li> ExpansionStrategy : The expansion strategy used to create/store new child(ren) from a
 *leaf node. Allows to decide whether a node should be expanded or not.</li>
 *<li> SimulationStrategy: The simulation strategy is used to select moves in self-play until the
 *end of the game. This is where you should put your simulation strategy.</li>
 * </ul>
 * @author lhm
 *
 * @note I could not manage to provide fully customizable mcts implementation and had to make
 *implementation choices :
 *<ul>
 *<li> Selection : Fixed using UCT (Upper Confidence Bound applied to Trees)</li>
 *<li> Expansion : Ideally, the user should be able to choose how he wants to expand a leaf node.
 *The current implementation is 'one node per simulated game + children of a node when its visit
 *count equals MCTS::_vis_expand_thresh (\see MCTS::set_expand_thresh()) </li>
 *<li> Final move selection : this implementation uses the "Robust child" (i.e. the most visited
 *one) but there are others (max child, robust-max child, secure child...) </li>
 */

// Standard headers
#include <algorithm>
#include <chrono>
#include <cmath>
#include <list>
#include <map>
#include <memory>
#include <ostream>
#include <random>
#include <set>
#include <type_traits>

#include <iostream>

#ifndef MCTS_HPP
#define MCTS_HPP

namespace mcts {

/*************************************************************************************************/
/*!
 * @brief Very basic stopwatch to measure ellapsed time since creation and tops.
 *
 * @note You do not have to use it directly. It could have lived under mcts::MCTS but it still might
 * be of help to have it...
 */
template<class Clock = std::chrono::high_resolution_clock, class Unit = std::chrono::milliseconds>
class stopwatch
{
public:
    explicit stopwatch() noexcept
      : _start{ Clock::now() }
      , _top{ _start }
    {}
    ~stopwatch() noexcept = default;

    Unit top() noexcept
    {
        auto now{ Clock::now() };
        Unit ret{ std::chrono::duration_cast<Unit>(now - _top) };
        _top = now;
        return ret;
    }

    Unit elapsed() const noexcept
    {
        return std::chrono::duration_cast<Unit>(Clock::now() - _start);
    }

private:
    typedef std::chrono::time_point<Clock> TimePt;
    const TimePt                           _start, _top;
};

/*************************************************************************************************/
/*!
 * @brief The Printable class provide an interface for printable objects.
 *
 * @note This is mainly used for debugging purposes.
 */
class Printable
{
protected:
    /*!
     * \brief print Override this to give your object a operator<< overloading (usefull when
     * debugging).
     * \param os An output stream object
     */
    virtual std::ostream& print(std::ostream& os) const noexcept { return os << this; }

    virtual ~Printable() noexcept = default;

private:
    friend std::ostream& operator<<(std::ostream& os, const Printable& obj) noexcept
    {
        return obj.print(os);
    }
};

/*************************************************************************************************/
/*!
 * @brief The State class represents a state of the game.
 * A state should be as light and easily-exploitable as possible.
 *
 * @note Derived classes should provide the following :
 * - operator<
 * - operator=
 * - copy constructor
 * @note An additional 'print' method can be implemented, mainly for debugging purposes.
 */
class State : public Printable
{
public:
    State() noexcept = default;
    virtual ~State() noexcept = default;
};

/*************************************************************************************************/
/*!
 * @brief The TerminalCriteria is the interface that allows to know whether a \a State is terminal
 * (i.e. no more \a Move can be produced from it).
 */
template<class St>
class TerminalCriteria
{
public:
    TerminalCriteria() noexcept = default;
    virtual ~TerminalCriteria() noexcept = default;

    /*!
     * \brief finished evaluates if the given state is terminal (i.e. the game is over)
     * \param state the state of the game
     * \return true if the state is terminal, false otherwise
     */
    virtual bool finished(const St& state) noexcept = 0;
};

/*************************************************************************************************/
/*!
 * @brief The TerminalEval is the interface that allows to evaluate a final \a State
 * @note Your evaluation should be positive when the final state is favorable to you, negative
 * otherwise.
 */
template<class St>
class TerminalEval
{
public:
    TerminalEval() noexcept = default;
    virtual ~TerminalEval() noexcept = default;

    /*!
     * \brief eval an evaluation function of a final state (\see TerminalCriteria) of the game.
     * \param state the state of the game
     * \return a value that should be positive if the result is good for the player, negative if it
     * is negative for the player, something else otherwise.
     */
    virtual float eval(const St& state) noexcept = 0;
};

/**
 * @brief Adapts the evaluation value during backp-ropagation
 */
template<class St>
class BackPropagationStrategy
{
public:
    explicit BackPropagationStrategy() noexcept = default;
    virtual ~BackPropagationStrategy() noexcept = default;

    /**
     * @param state the state we adjust the evaluation value for
     * @param eval The evaluation value backpropagated (\see TerminalEval)
     * @return The modified evaluation value for the current state
     */
    virtual float adjust(const St& state, float eval) = 0;
};

/*************************************************************************************************/
/*!
 * @brief The Move class represents a move of the game.
 * It can be applied to a \a State in order to modify it.
 * @note Derived classes should be default constructible and have an operator== implementation
 */
template<class St>
class Move : public Printable
{
public:
    Move() noexcept = default;
    virtual ~Move() noexcept = default;

    /*!
     * @brief apply Apply the move to a given \a State.
     * @param state the state to apply the move to.
     * @return true if the application of the move changed the state of the game, false otherwise.
     * @note This should modify the state of the game, unless the move is forbidden for the given
     * state.
     */
    [[maybe_unused]] virtual bool apply(St& state) const noexcept = 0;
};

/*************************************************************************************************/
/*!
 * @brief The expansion strategy controls how leaf nodes/states are expended.
 */
template<class St, class Mv>
class ExpansionStrategy
{
public:
    explicit ExpansionStrategy() noexcept = default;
    virtual ~ExpansionStrategy() noexcept = default;

    /*!
     * \brief expand Produces one \a Move from a given \a State.
     * The move will lead to the creation/storage of a new node.
     */
    virtual Mv expand(const St& state) noexcept = 0;

    /*!
     * \brief expandAll produces every possible move from the current state
     */
    virtual std::set<Mv> expandAll(const St& state) noexcept = 0;

    /*!
     * \brief is_expandable evaluates if the given state is expandable
     * \param state the state of the game
     * \return true if the state is expandable (i.e. there are moves left).
     */
    virtual bool is_expandable(const St& state) noexcept = 0;
};

/*************************************************************************************************/
/*!
 * @brief
 */
template<class St, class Mv>
class SimulationStrategy
{
public:
    explicit SimulationStrategy() noexcept = default;
    virtual ~SimulationStrategy() noexcept = default;

    /*!
     * \brief simulate creates a \a Move that will be performed on the current state to
     * simulate a game.
     * \return the move produced by the simulation.
     */
    virtual Mv simulate(const St& state) noexcept = 0;
};

/*************************************************************************************************/
/*!
 * @brief The Node class represents a node of the Monte-Carlo Tree.
 */
template<class St, class Mv>
class Node
{
public:
    using N = Node<St, Mv>;

public:
    /*!
     * @brief Node
     * @param state The state associated to the node
     * @param mv The move that led to that node
     * @param term The termination criteria (used to know it the node is final)
     * @param exp The expansion strategy (used to know potential children)
     * @param parent The parent of the node
     */
    Node(const St&                                         state,
         const Mv&                                         mv,
         const std::shared_ptr<TerminalCriteria<St>>&      term,
         const std::shared_ptr<ExpansionStrategy<St, Mv>>& exp,
         N*                                                parent = nullptr) noexcept
      : _state{ state }
      , _move{ mv }
      , _parent{ parent }
      , _terminal{ term->finished(_state) }
      , _unexplored{}
    {
        if (!_terminal)
            _unexplored = exp->expandAll(_state);

        if (nullptr != _parent)
            _parent->add_child(this);
    }
    ~Node() noexcept
    {
        for (auto& nxt : _children)
            delete (nxt);
    }

    Node(const Node&) noexcept = delete;
    Node operator=(const Node&) noexcept = delete;

    void set_parent(N* node) noexcept { _parent = node; }
    bool is_terminal(void) const noexcept { return _terminal; }

    bool has_children(void) const noexcept { return !std::empty(_children); }

    /*!
     * \brief is_expandable evaluates if the node is expandable
     * \return true if it is a expandable, false otherwise.
     */
    bool is_expandable(void) const noexcept { return !std::empty(_unexplored); }

    /*!
     * \brief add_child add a child to the node
     * \param child the child to be added
     */
    void add_child(N* child) noexcept
    {
        if (_unexplored.count(child->move())) {
            _unexplored.erase(child->move());
            _children.push_back(child);
        }
    }

    /*!
     * \brief grab_child get a child  of the node. It will give you the ownership over that child
     * (and the responsability to delete it). After that call, the corresponding node is not a child
     * of the current node anymore.
     * \param mv the move
     * \return the child node if it exists, nullptr otherwise.
     */
    N* grab_child(const Mv& mv) noexcept
    {
        N* ret{ nullptr };
        for (auto& child : _children) {
            if (mv == child->move()) {
                ret = child;
                _children.remove(child);
                break;
            }
        }
        if (nullptr != ret) {
            ret->set_parent(nullptr);
        }
        return ret;
    }

    N* random_child() noexcept
    {
        if (std::empty(_children))
            return nullptr;

        std::vector<N*> ret{};
        std::sample(std::begin(_children),
                    std::end(_children),
                    std::back_inserter(ret),
                    1,
                    std::mt19937{ std::random_device{}() });

        return ret[0];
    }

    /*!
     * \brief random_move will randomly get a move from the possible unexplored moves.
     * \return a Move
     * \note make sure the Node is a leaf before calling this method
     */
    Mv random_move() noexcept
    {
        if (std::empty(_unexplored))
            return {};

        std::vector<Mv> ret{};
        std::sample(std::begin(_unexplored),
                    std::end(_unexplored),
                    std::back_inserter(ret),
                    1,
                    std::mt19937{ std::random_device{}() });

        return ret[0];
    }

    /*!
     * \brief update update the statistics associated to the node with the given result.
     * \param res the result (obtained from a simulation).
     */
    void update(float res) noexcept
    {
        ++_visit_count;
        _res_count += res;
        _val = _res_count / _visit_count;
    }

    /**
     * Getters
     */

    const St&   state(void) noexcept { return _state; }
    const Mv&   move(void) noexcept { return _move; }
    auto        parent(void) noexcept { return _parent; }
    const auto& children(void) noexcept { return _children; }
    const auto& unexploredMoves(void) noexcept { return _unexplored; }

    const auto& visit_count(void) noexcept { return _visit_count; }
    const auto& res_count(void) noexcept { return _res_count; }
    const auto& val(void) noexcept { return _val; }

private:
    /**
     * Tree related parameters
     */
    const St      _state;      /*< The game state associated to the node */
    const Mv      _move;       /*< The move that led to that state */
    N*            _parent;     /*< The node that produced it */
    bool          _terminal;   /*< The node corresponds to a terminal state */
    std::set<Mv>  _unexplored; /*< The move that have not been explored yet */
    std::list<N*> _children;   /*< Children emplaced (either by simulation or expansion) */

    /**
     * Stats related parameters
     * updated during back-propagation
     */
    uint32_t _visit_count{ 0 }; /*< The number of times the node has been visited */
    float    _res_count{ 0 };   /*< The total of every results backpropagated to this node */
    float    _val{ 0 };         /*< The value computed for the node - updated by backPropagation */
};

template<class St, class Mv>
class MCTS
{
public:
    using N = Node<St, Mv>;

public:
    MCTS(const St&                                           initialState,
         const std::shared_ptr<SimulationStrategy<St, Mv>>&  simulationStrategy,
         const std::shared_ptr<ExpansionStrategy<St, Mv>>&   expansionStrategy,
         const std::shared_ptr<TerminalCriteria<St>>&        terminalCriteria,
         const std::shared_ptr<TerminalEval<St>>&            terminalEval,
         const std::shared_ptr<BackPropagationStrategy<St>>& backpropStrategy)
    noexcept
      : _sim{ simulationStrategy }
      , _exp{ expansionStrategy }
      , _termCrit{ terminalCriteria }
      , _termEval{ terminalEval }
      , _backProp{ backpropStrategy }
      , _root{ new N(initialState, {}, _termCrit, _exp) }
    {}

    ~MCTS() noexcept
    {
        if (nullptr != _root)
            delete (_root);
    }

    /*!
     * \brief advance updates the tree with the given move
     * This will free previously allocated parts of the tree that are no longer relevant.
     * \param mv the move to update the tree with
     * \return true in case of success, false otherwise
     */
    [[maybe_unused]] bool advance(const Mv& mv) noexcept
    {
        if (nullptr == _root)
            return false;

        auto next_root{ _root->grab_child(mv) };
        if (nullptr == next_root) {
            auto st{ _root->state() };
            mv.apply(st);
            next_root = new N(st, mv, _termCrit, _exp);
        }

        delete (_root);
        _root = next_root;

        return true;
    }

    /*!
     * @brief compute Perform the MCTS algorithm
     * \return The next move to play
     */
    Mv compute(void) noexcept
    {
        if (nullptr == _root || nullptr == _sim || nullptr == _exp || nullptr == _termCrit ||
            nullptr == _termEval)
            return {};

        /**
         * while(there is time left) {
         *  4 steps (selection, expansion, simulation, backpropagation)
         * }
         *
         * best Mv = argmax( n in root's children nodes)(N.visit_count)
         */
        stopwatch sw;
        // uint32_t  cur_it{ 0 };
        while (_time_max > sw.elapsed() /*&& (++cur_it < _max_it)*/) {
            N* cur_node{ _root };

            /**
             * Selection
             * Recursively select nodes until a leaf is found
             *
             * If visit count > _vis_uct_thresh
             *     Select the node n in (reachable from cur_node)
             *     that maximizes UCT.
             * else use random selection
             */
            cur_node = _select(cur_node);
            if (cur_node->is_terminal()) {
                _backpropagate(cur_node, cur_node->state());
                continue;
            }

            /**
             * Expansion
             *
             * - Expand the first node that is not in the tree
             * - Also expand all the children of a node when its visit_count == _vis_expand_thresh
             */
            cur_node = _expand(cur_node);

            /**
             * Simulation
             */
            auto st{ _simulate(cur_node) };

            /**
             * BackPropagation
             *
             * Propagate the result of the simulation backwards from the leaf node to the root.
             */
            _backpropagate(cur_node, st);
        }

        /**
         * Final move selection
         *
         * The best move is the root child that maximizes the visit_count
         */
        return _bestMoveSelection();
    }

    /**
     * Setters
     */
    void set_time_max(uint32_t t) noexcept { _time_max = std::chrono::milliseconds(t); }
    void set_default_c(float c) noexcept { _C = c; }
    void set_expand_thresh(uint32_t thresh) noexcept { _vis_expand_thresh = thresh; }
    void set_uct_thresh(uint32_t thresh) noexcept { _vis_uct_thresh = thresh; }
    void set_max_it(uint32_t max_it) noexcept { _max_it = max_it; }

private:
    /**
     * 4 stages implementations
     * + Final move selection
     */

    /*!
     * \brief _select implementation of the "selection stage"
     * \param n the initial node, will be modified during execution to point to the selcted node
     */
    N* _select(N* n) noexcept
    {
        N* ret{ n };
        while (!ret->is_expandable() && ret->has_children()) {
            if (ret->visit_count() > _vis_uct_thresh) {
                N*    sel_node{ nullptr };
                float max_uct{ -std::numeric_limits<float>::max() };
                for (const auto& nxt : ret->children()) {
                    float uct{ ret->val() +
                               _C * (float)sqrt(log(ret->visit_count()) / nxt->visit_count()) };
                    if (uct > max_uct) {
                        sel_node = nxt;
                        max_uct = uct;
                    }
                }
                ret = sel_node;
            } else {
                ret = ret->random_child();
            }
        }
        return ret;
    }

    /*!
     * \brief _expand implementation of the "expansion stage"
     * \param n the initial node to be expanded, its children will be modified
     * \return a new node if expansion took place, nullptr otherwise
     */
    N* _expand(N* n) noexcept
    {
        N* ret{ n };

        if (!ret->is_expandable())
            return ret;

        if (ret->visit_count() >= _vis_expand_thresh) {
            auto new_st{ n->state() };
            auto new_mv{ n->random_move() };
            new_mv.apply(new_st);
            ret = new N(new_st, new_mv, _termCrit, _exp, n);
        }

        /*
        if (ret->visit_count() >= _vis_expand_thresh) {
            // expand all possible moves
            auto unexploredMoves{ ret->unexploredMoves() };
            for (const auto& new_mv : unexploredMoves) {
                auto new_st{ ret->state() };
                new_mv.apply(new_st);
                ret = new N(new_st, new_mv, _termCrit, _exp, n);
            }
        } else {
            // expand a random move
            auto new_st{ n->state() };
            auto new_mv{ n->random_move() };
            new_mv.apply(new_st);
            ret = new N(new_st, new_mv, _termCrit, _exp, n);
        }
        */

        return ret;
    }

    /*!
     * \brief _simulate implementation of the "simulation stage"
     * \param n the starting node of the simulation
     * \return the final state of the simulation
     */
    St _simulate(N* n) noexcept
    {
        St ret{ n->state() };
        Mv mv{};
        while (!_termCrit->finished(ret)) {
            mv = _sim->simulate(ret);
            mv.apply(ret);
        }

        return ret;
    }

    /*!
     * \brief _backpropagate implementation of the "backpropagation stage"
     * \param n the initial leaf node to start propagating from, will be modified during execution
     * to point to the initial (root) node.
     * \param state the final state produced by the simulation (or the selection).
     */
    void _backpropagate(N* n, const St& state) noexcept
    {
        auto eval{ _termEval->eval(state) };
        auto curNode{ n->parent() };

        n->update(_backProp->adjust(n->state(), eval));

        while (nullptr != curNode) {
            curNode->update(_backProp->adjust(curNode->state(), eval));
            curNode = curNode->parent();
        }

        n = curNode;
    }

    /*!
     * \brief _bestMoveSelection Selects the best move to play (best child of the root node)
     * \return the best move
     */
    Mv _bestMoveSelection(void) noexcept
    {
        // Debug
        // N* retNode{ nullptr };

        Mv    ret;
        float max_val{ -std::numeric_limits<float>::max() };
        for (const auto& n : _root->children()) {
            if (n->val() >= max_val) {
                max_val = n->val();
                ret = n->move();
            }
        }

        // TODO - If no candidate, play random move ?

        /*
        std::cout << "Best move : (vis, wins, val): " << retNode->visit_count() << ", "
                  << retNode->res_count() << ", " << retNode->val() << std::endl;
                  */

        return ret;
    }

protected:
    /**
     * User provided interfaces
     */

    std::shared_ptr<SimulationStrategy<St, Mv>>  _sim;
    std::shared_ptr<ExpansionStrategy<St, Mv>>   _exp;
    std::shared_ptr<TerminalCriteria<St>>        _termCrit;
    std::shared_ptr<TerminalEval<St>>            _termEval;
    std::shared_ptr<BackPropagationStrategy<St>> _backProp;

private:
    /**
     * Tree related structures
     */
    N* _root;

    /**
     * Customizable params
     */

    std::chrono::milliseconds _time_max{ default_time };
    float                     _C{ default_c };
    uint32_t                  _vis_expand_thresh{ default_vis };
    uint32_t                  _vis_uct_thresh{ default_vis_thresh };
    uint32_t                  _max_it{ default_max_it };
    /**
     * Other params
     */

private:
    /**
     * Default values for customizable params
     */

    static constexpr uint32_t default_time{ 1000 };      /*< time in ms */
    static constexpr float    default_c{ 0.7 };          /*< UCT constant */
    static constexpr uint32_t default_vis{ 5 };          /*< nb of visits before expansion */
    static constexpr uint32_t default_vis_thresh{ 5 };   /*< Do not apply UCT if vis_count < this */
    static constexpr uint32_t default_max_it{ 1000000 }; /*< Do not make more than this nb of it */
};

} // namespace mcts

#endif // MCTS_HPP