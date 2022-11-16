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
#include <chrono>
#include <cmath>
#include <list>
#include <map>
#include <memory>
#include <ostream>
#include <type_traits>

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
    virtual int32_t eval(const St& state) noexcept = 0;
};

/*************************************************************************************************/
/*!
 * @brief The Move class represents a move of the game.
 * It can be applied to a \a State in order to modify it.
 * @note Derived classes should be default constructible.
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
    [[maybe_unused]] virtual bool apply(St& state) noexcept = 0;
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
     * \note The produced move must not lead to already expanded nodes.
     * To ensure this, you should make sure that subsequent calls to expand (on the same state) will
     * lead to different moves.
     */
    virtual Mv expand(const St& state) noexcept = 0;

    /*!
     * \brief is_expandable evaluates if the given state is expandable
     * \param state the state of the game
     * \return true if the state is expandable (i.e. it is not final and there are moves left that
     * the \a ExpansionStrategy never returned).
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
     * @param parent The parent of the node
     */
    Node(const St& state, const Mv& mv, const N* parent = nullptr) noexcept
      : _state{ state }
      , _move{ mv }
      , _parent{ parent }
    {
        if (nullptr != _parent)
            _parent->add_child(this);
    }
    ~Node() noexcept
    {
        if (!is_leaf())
            for (auto& nxt : _children)
                delete (nxt);
    }

    /*!
     * \brief is_leaf evaluates if the node is a leaf
     * \return true if it is a leaf, false otherwise.
     */
    bool is_leaf(void) const noexcept { return std::empty(_children); }

    /*!
     * \brief add_child add a child to the node
     * \param child the child to be added
     */
    void add_child(N* child) noexcept { _children.push_back(child); }

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
        return ret;
    }

    /*!
     * \brief update update the statistics associated to the node with the given result.
     * \param res the result (obtained from a simulation).
     */
    void update(int32_t res) noexcept
    {
        ++_visit_count;
        _res_count += res;
        _val = (float)_res_count / _visit_count;
    }

    /**
     * Getters
     */

    const auto& state(void) noexcept { return _state; }
    const auto& move(void) noexcept { return _move; }
    auto        parent(void) noexcept { return _parent; }
    const auto& children(void) noexcept { return _children; }

    const auto& visit_count(void) noexcept { return _visit_count; }
    const auto& res_count(void) noexcept { return _res_count; }
    const auto& val(void) noexcept { return _val; }

private:
    /**
     * Tree related parameters
     */
    const St      _state;    /*< The game state associated to the node */
    const Mv      _move;     /*< The move that led to that state */
    const N*      _parent;   /*< The node that produced it */
    std::list<N*> _children; /*< Children emplaced (either by simulation or expansion) */

    /**
     * Stats related parameters
     * updated during back-propagation
     */
    uint32_t _visit_count{ 0 }; /*< The number of times the node has been visited */
    uint32_t _res_count{ 0 };   /*< The total of every results backpropagated to this node */
    float    _val{ 0 };         /*< The value computed for the node - updated by backPropagation */
};

template<class St, class Mv>
class MCTS
{
public:
    using N = Node<St, Mv>;

public:
    MCTS(const St&                                          initialState,
         const std::shared_ptr<SimulationStrategy<St, Mv>>& simulationStrategy,
         const std::shared_ptr<ExpansionStrategy<St, Mv>>&  expansionStrategy,
         const std::shared_ptr<TerminalCriteria<St>>&       terminalCriteria,
         const std::shared_ptr<TerminalEval<St>>&           terminalEval)
    noexcept
      : _root{ new N(initialState, Mv()) }
      , _sim{ simulationStrategy }
      , _exp{ expansionStrategy }
      , _termCrit{ terminalCriteria }
      , _termEval{ terminalEval }
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
        if (nullptr == next_root)
            return false;

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
        while (_time_max > sw.elapsed()) {
            N* cur_node{ _root };

            /**
             * Selection
             * Recursively select nodes until a leaf is found
             *
             * If visit count > _vis_uct_thresh
             *     Select the node n in (reachable from cur_node)
             *     that maximizes UCT.
             * else use the simulation strategy (not impl)
             */
            _select(cur_node);

            /**
             * Expansion
             *
             * - Expand the first node that is not in the tree
             * - Also expand all the children of a node when its visit_count == _vis_expand_thresh
             */
            _expand(cur_node);

            /**
             * Simulation
             */
            _simulated(cur_node);

            /**
             * BackPropagation
             *
             * Propagate the result of the simulation backwards from the leaf node to the root.
             */
            _backpropagate(cur_node);
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

private:
    /**
     * 4 stages implementations
     * + Final move selection
     */

    /*!
     * \brief _select implementation of the "selection stage"
     * \param n the initial node, will be modified during execution to point to the selcted node
     */
    void _select(N* n) noexcept
    {
        if (n->visit_count() > _vis_uct_thresh) {
            while (!n->is_leaf()) {
                float max_uct{ 0 };
                N*    sel_node{ nullptr };
                for (const auto& nxt : n->children()) {
                    float uct{ n->val() +
                               _C * (float)sqrt(log(n->visit_count()) / nxt->visit_count()) };
                    if (uct >= max_uct) {
                        sel_node = nxt;
                        max_uct = uct;
                    }
                }
                n = sel_node;
            }
        }
    }

    /*!
     * \brief _expand implementation of the "expansion stage"
     * \param n the initial node, will be modified during execution to point to the expanded node
     */
    void _expand(N* n) noexcept
    {
        const auto& st{ St(n->state()) };
        auto        mv{ Mv() };
        if (n->visit_count() == _vis_expand_thresh) {
            // expand all
            while (_exp->is_expandable(st)) {
                auto new_st{ st };
                mv = _exp->expand(new_st);
                mv.apply(new_st);
                n = new N(new_st, mv, n);
            }
        } else {
            // expand the first
            if (_exp->is_expandable(st)) {
                auto new_st{ st };
                mv = _exp->expand(new_st);
                mv.apply(new_st);
                n = new N(new_st, mv, n);
            }
        }
    }

    /*!
     * \brief _simulate implementation of the "simulation stage"
     * \param n the initial node, will be modified during execution to point to the node of the
     * final simulation move.
     */
    void _simulate(N* n) noexcept
    {
        auto st{ St(n->state()) };
        auto mv{ Mv() };
        while (!_termCrit->finished(st)) {
            // Create next node from the move
            mv = _sim->simulate(st);
            mv.apply(st);
            n = new N(st, mv, n);
        }
    }

    /*!
     * \brief _backpropagate implementation of the "backpropagation stage"
     * \param n the initial leaf node to start propagating from, will be modified during execution
     * to point to the initial (root) node.
     */
    void _backpropagate(N* n) noexcept
    {
        auto eval{ _termEval->eval(n->state()) };
        while (nullptr != n->parent()) {
            n->update(eval);
            n = n->parent();
        }
    }

    /*!
     * \brief _bestMoveSelection Selects the best move to play (best child of the root node)
     * \return the best move
     */
    Mv _bestMoveSelection(void) noexcept
    {
        Mv       ret;
        uint32_t max_visit_count{ 0 };
        for (const auto& n : _root->children()) {
            if (n->visit_count() >= max_visit_count) {
                max_visit_count = n->visit_count();
                ret = n->move();
            }
        }

        return ret;
    }

protected:
    /**
     * User provided interfaces
     */

    std::shared_ptr<SimulationStrategy<St, Mv>> _sim;
    std::shared_ptr<ExpansionStrategy<St, Mv>>  _exp;
    std::shared_ptr<TerminalCriteria<St>>       _termCrit;
    std::shared_ptr<TerminalEval<St>>           _termEval;

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

    /**
     * Other params
     */

private:
    /**
     * Default values for customizable params
     */

    static constexpr uint32_t default_time{ 1000 };     /*< time in ms */
    static constexpr float    default_c{ 0.5 };         /*< UCT constant */
    static constexpr uint32_t default_vis{ 7 };         /*< nb of visits before expansion */
    static constexpr uint32_t default_vis_thresh{ 30 }; /*< Do not apply UCT if vis_count < this */
};

} // namespace mcts

#endif // MCTS_HPP
