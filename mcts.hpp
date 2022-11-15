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
 *<li> ExpansionStrategy : The expansion strategy used to create/store new child(ren) from a leaf
 *node. Allows to decide whether a node should be expanded or not.</li>
 *<li> SimulationStrategy : The simulation strategy is used to select moves in self-play until the
 *end of the game. This is where you should put your simulation strategy.</li>
 * </ul>
 * @author lhm
 *
 * @note I could not manage to provide fully customizable mcts implementation and had to make
 *implementation choices :
 * <ul>
 * <li> Selection : Fixed using UCT (Upper Confidence Bound applied to Trees)</li>
 * <li> Expansion : Ideally, the user should be able to choose how he wants to expand a leaf node.
 * (One node ? more ? all possible children ?) and know whether its expansion is not alreay stored
 *in the tree. </li>
 * <li> Final move selection : this implementation uses the "Robust child", but there are others
 *(max child, robust-max child, secure child...) </li>
 */

// Standard headers
#include <chrono>
#include <cmath>
#include <map>
#include <ostream>
#include <type_traits>

#ifndef MCTS_HPP
#define MCTS_HPP

namespace mcts {

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
    virtual ~State() noexcept = default;

    /*!
     * \brief eval an evaluation function of the current state of the game.
     * \return a value that should be 0 if the state is not final,
     * 1 if the result is good for the player, -1 if it is negative for the player, 0.5 otherwise.
     */
    virtual int32_t eval(void) noexcept = 0;
};

/*************************************************************************************************/
/*!
 * @brief The Move class represents a move of the game.
 * It can be applied to a \a State in order to modify it.
 * @note Derived classes should be default constructible.
 */
template<class St, typename = std::enable_if_t<std::is_base_ofase_of_v<State, St>>>
class Move : public Printable
{
public:
    operator bool() noexcept { return _possible; }

public:
    Move(bool possible = false) noexcept
      : _possible{ possible }
    {}
    virtual ~Move() noexcept = default;

    /*!
     * @brief apply Apply the move to a given \a State.
     * @param state the state to apply the move to.
     *
     * @note This should modify the state of the game, unless the move is forbidden for the given
     * state.
     */
    virtual void apply(St& state) noexcept = 0;

protected:
    bool _possible;
};

/*************************************************************************************************/
struct Stat
{
    uint32_t visit_count{ 0 }; /*< The number of times the node has been visited */
    uint32_t res_count{ 0 };   /*< The total of every results backpropagated to this node */
    float    val{ 0 };         /*< The value computed for the node - updated by backPropagation */

    void update(int32_t res) noexcept
    {
        ++visit_count;
        res_count += res;
        val = (float)res_count / visit_count;
    }
};

/*************************************************************************************************/
/*!
 * @brief The expansion strategy controls how leaf nodes/states are expended.
 */
template<class St,
         class Mv,
         typename = std::enable_if_t<std::is_base_of_v<State, St>>,
         typename = std::enable_if_t<std::is_base_of_v<Move, Mv>>>
class ExpansionStrategy
{
public:
    /*!
     * \brief expand Produces one \a Move from a given \a State.
     * The move will lead to the creation/storage of a new node.
     * \note The produced move must not lead to already expanded nodes.
     * To ensure this, you should make sure that subsequent calls to expand (on the same state) will
     * lead to different moves.
     *
     * In case there is nothing to expand (i.e. final State), you should create an impossible move.
     * (\see Move::_possible)
     */
    virtual Mv expand(const St& state) noexcept = 0;

    explicit ExpansionStrategy() noexcept = default;
    virtual ~ExpansionStrategy() noexcept = default;
};

/*************************************************************************************************/
/*!
 * @brief TODO
 */
template<class St,
         class Mv,
         typename = std::enable_if_t<std::is_base_of_v<State, St>>,
         typename = std::enable_if_t<std::is_base_of_v<Move, Mv>>>
class SimulationStrategy
{
public:
    /*!
     * \brief simulate creates a \a Move that will be performed on the current state to
     * simulate a game.
     * \return the move produced by the simulation. It there are no move to perform (i.e. the state
     * is final), it should return an impossible move (\see Move::_possible)
     */
    virtual Mv simulate(const St& state) noexcept = 0;

    explicit SimulationStrategy() noexcept = default;
    virtual ~SimulationStrategy() noexcept = default;
};

/*************************************************************************************************/
/*!
 * @brief The Node class represents a node of the Monte-Carlo Tree.
 */
template<class St,
         class Mv,
         typename = std::enable_if_t<std::is_base_of_v<State, St>>,
         typename = std::enable_if_t<std::is_base_of_v<Move, Mv>>>
class Node
{
public:
    using N = Node<St, Mv>;

public:
    /*!
     * @brief Node
     * @param state The state associated to the node
     * @param move The move that led to that node
     * @param parent The parent of the node
     */
    Node(const St& state, const Mv& move, const N* parent = nullptr) noexcept
      : _state{ state }
      , _move{ move }
      , _parent{ parent }
      , _stats{}
    {
        if (nullptr != _parent)
            _parent->add_child(this);
    }
    ~Node() noexcept = default;

    const auto& state(void) noexcept { return _state; }
    const auto& move(void) noexcept { return _move; }
    const auto  parent(void) noexcept { return _parent; }
    const auto& children(void) noexcept { return _children; }
    auto&       stats(void) noexcept { return _stats; }

    void add_child(N* child) noexcept { _children[child->move()] = child; }

protected:
private:
    const St         _state;    /*< The game state associated to the node */
    const Mv         _move;     /*< The move that led to that state */
    const N*         _parent;   /*< The node that produced it */
    std::map<Mv, N*> _children; /*< Children emplaced (either by simulation or expansion) */

    Stat _stats; /*< The statistics associated to this node (updated by every step) */
};

template<class St,
         class Mv,
         typename = std::enable_if_t<std::is_base_of_v<State, St>>,
         typename = std::enable_if_t<std::is_base_of_v<Move, Mv>>>
class MCTS
{
public:
    MCTS(const St&                     initialState,
         const SimulationStrategy<St>& simulationStrategy,
         const ExpansionStrategy<St>&  expansionStrategy)
    noexcept
      : _root{ N(initialState, Mv()) }
      , _sim{ simulationStrategy }
      , _exp{ expansionStrategy }
    {
        _nodes[initialState] = &_root;
    }

    ~MCTS() noexcept
    {
        for (auto& [st, n] : nodes) {
            if (nullptr != n) {
                delete (n);
                n = nullptr;
            }
        }
        _nodes.clear();
    }

    /*!
     * @brief compute Perform the MCTS algorithm
     * \return The next move to play
     */
    Mv compute(void) noexcept
    {
        /**
         * while(time) {
         *  4 steps (selection, expansion, simulation, backpropagation)
         * }
         *
         * best Mv = argmax( n in root's children nodes)(N.visit_count)
         */
        stopwatch sw;
        while (_time_max > sw.elapsed()) {
            N* cur_node{ &_root };

            /**
             * Selection
             * Recursively select nodes until a leaf is found
             *
             * If visit count > _vis_uct_thresh
             *     Select the node n in (reachable from cur_node)
             *     that maximizes UCT.
             * else use the simulation strategy
             */
            {
                if (cur_node->stats().visit_count > _vis_uct_thresh) {
                    while (!std::empty(cur_node->children())) {
                        float max_uct{ 0 };
                        N*    sel_node{ nullptr };
                        for (const auto& [mv, n] : cur_node->children()) {
                            float uct{ cur_node->stats().val +
                                       _C * (float)sqrt(log(cur_node->stats().visit_count) /
                                                        n->stats().visit_count) };
                            if (uct >= max_uct) {
                                sel_node = n;
                                max_uct = uct;
                            }
                        }
                        cur_node = sel_node;
                    }
                }
            }

            /**
             * Expansion
             *
             * - Expand the first node that is not in the tree
             * - Also expand all the children of a node when its visit_count == _vis_expand_thresh
             */
            {
                auto st{ St(cur_node->state()) };
                auto mv{ Mv() };
                if (cur_node->stats().visit_count == _vis_expand_thresh) {
                    // expand all
                    while ((mv = _exp.expand(st))) {
                        auto new_st{ st };
                        mv.apply(new_st);
                        _nodes[new_st] = new N(new_st, mv, cur_node);
                        cur_node = _nodes[new_st];
                    }
                } else {
                    // expand the first
                    mv = _exp.expand(st);
                    if (mv) {
                        auto new_st{ st };
                        mv.apply(new_st);
                        _nodes[new_st] = new N(new_st, mv, cur_node);
                        cur_node = _nodes[new_st];
                    }
                }
            }

            /**
             * Simulation
             */
            {
                auto st{ St(cur_node->state()) };
                auto mv{ Mv() };
                while ((mv = _sim.simulate(cur_node->state()))) {
                    // Create next node from the move
                    mv.apply(st);
                    _nodes[st] = new N(st, mv, cur_node);
                    cur_node = _nodes[st];
                }
            }

            /**
             * BackPropagation
             *
             * Propagate the result of the simulation backwards from the leaf node to the root.
             */
            auto eval{ st.eval() };
            while (nullptr != cur_node) {
                cur_node->stats().update(eval);
                cur_node = cur_node->parent();
            }
        }

        /**
         * Final move selection
         *
         * The best move is the root child that maximizes the visit_count
         */
        {
            Mv       ret;
            uint32_t max_visit_count{ 0 };
            for (const auto& [mv, n] : _root.children()) {
                if (n->stats().visit_count >= max_visit_count) {
                    max_visit_count = n->stats().visit_count;
                    ret = mv;
                }
            }

            return ret;
        }
    }

    // Setters
    void     set_time_max(uint32_t t) noexcept { _time_max = std::chrono::milliseconds(t); }
    void     set_default_c(float c) noexcept { _C = c; }
    uint32_t set_expand_thresh(uint32_t thresh) noexcept { _vis_expand_thresh = thresh; }
    uint32_t set_uct_thresh(uint32_t thresh) noexcept { _vis_uct_thresh = thresh; }

protected:
    SimulationStrategy<St> _sim;
    ExpansionStrategy<St>  _exp;

private:
    // Tree related structures
    N                _root;  /*< The root of the tree */
    std::map<St, N*> _nodes; /*< Container for every nodes of the tree */

    // Customizable params
    std::chrono::milliseconds _time_max{ default_time };
    float                     _C{ default_c };
    uint32_t                  _vis_expand_thresh{ default_vis };
    uint32_t                  _vis_uct_thresh{ default_vis_thresh };

    // Other params

private:
    static constexpr uint32_t default_time{ 1000 };     /*< time in ms */
    static constexpr float    default_c{ 0.5 };         /*< UCT constant */
    static constexpr uint32_t default_vis{ 7 };         /*< nb of visits before expansion */
    static constexpr uint32_t default_vis_thresh{ 30 }; /*< Do not apply UCT if vis_count < this */

private:
    /*!
     * @brief Very basic stopwatch to measure ellapsed time since creation and tops.
     */
    template<class Clock = std::chrono::high_resolution_clock,
             class Unit = std::chrono::milliseconds>
    class stopwatch
    {
    public:
        explicit stopwatch() noexcept
          : _start{ Clock::now() }
          , _top{ _start }
        {}

        Unit top() noexcept
        {
            auto now{ Clock::now() };
            Unit ret{ std::chrono::duration_cast<Unit>(now - top) };
            top = now;
            return ret;
        }

        Unit elapsed() const noexcept
        {
            return std::chrono::duration_cast<Unit>(Clock::now() - _start);
        }

    private:
        typedef std::chrono::time_point<Clock> TimePt;
        const TimePt                           _start, top;
    };
};

} // namespace mcts

#endif // MCTS_HPP
