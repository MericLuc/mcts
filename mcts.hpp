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
 *<li> Strategies : The strategies to perform the 4 strategic tasks of the MCTS algorithm.
 *<ul> <li> SelectionStrategy : The selection strategy that is recursively applied
 *untill a leaf node is found. Allows to select one of the children of a given node.</li>
 *<li> ExpansionStrategy : The expansion strategy used to create/store new child(ren) from a leaf
 *node. Allows to decide whether a node should be expanded or not.</li>
 *<li> SimulationStrategy : The simulation strategy is used to select moves in self-play untill the
 *end of the game. This is where you should put your simulation strategy.</li>
 *<li> BackPropagationStrategy : The backpropagation strategy propagates the result of the
 *simulation from leaf nodes to parent nodes.</li>
 * </ul>
 * </li>
 *</ul>
 * @author lhm
 */

// Standard headers
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
 * @note An additional 'print' method can be implemented, mainly for debugging purposes.
 */
class State : public Printable
{
public:
    virtual ~State() noexcept = default;

    /*!
     * \brief eval an evaluation function of the current state of the game.
     * \return a value that should be 0 if the state is not final (or if the game is null),
     * positive if the result is good for the player, negative otherwise.
     */
    virtual float eval(void) noexcept = 0;
};

/*************************************************************************************************/
/*!
 * @brief The Move class represents a move of the game.
 * It can be applied to a \a State in order to modify it.
 */
template<class St, typename = std::enable_if_t<std::is_base_of_v<State, St>>>
class Move : public Printable
{
public:
    operator bool() noexcept { return _possible; }

public:
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
/*!
 * @brief The Strategy class is the interface for every MCTS strategies.
 *
 * Basically, any strategy generate \a Action from a given \a State.
 */
template<class St, typename = std::enable_if_t<std::is_base_of_v<State, St>>>
class Strategy
{
public:
    explicit Strategy(St* state) noexcept
      : _state{ state }
    {}
    virtual ~Strategy() noexcept = default;

protected:
    St* _state;
};

/*************************************************************************************************/
/*!
 * @brief TODO
 */
template<class St,
         class Mv,
         typename = std::enable_if_t<std::is_base_of_v<State, St>>,
         typename = std::enable_if_t<std::is_base_of_v<Move, Mv>>>
class SelectionStrategy : public Strategy<St>
{
public:
    /*!
     * @brief select Produces a \a Move from the current \a State (Strategy::_state).
     * The move will lead to the selection of an existing node.
     */
    virtual void select(Mv& move) noexcept = 0;

    explicit SelectionStrategy(St* state) noexcept
      : Strategy<St>(state)
    {}
    virtual ~SelectionStrategy() noexcept = default;

protected:
private:
};

/*************************************************************************************************/
/*!
 * @brief TODO
 */
template<class St,
         class Mv,
         typename = std::enable_if_t<std::is_base_of_v<State, St>>,
         typename = std::enable_if_t<std::is_base_of_v<Move, Mv>>>
class ExpansionStrategy : public Strategy<St>
{
public:
    /*!
     * \brief expand Produces a \a Move from the current \a State (Strategy::_state).
     * The move will lead to the creation/storage of a new node.
     * \param move output move. Cast to bool to see if this is a possible move.
     */
    virtual void expand(Mv& move) noexcept = 0;

    explicit ExpansionStrategy(St* state) noexcept
      : Strategy<St>(state)
    {}
    virtual ~ExpansionStrategy() noexcept = default;

protected:
private:
};

/*************************************************************************************************/
/*!
 * @brief TODO
 */
template<class St,
         class Mv,
         typename = std::enable_if_t<std::is_base_of_v<State, St>>,
         typename = std::enable_if_t<std::is_base_of_v<Move, Mv>>>
class SimulationStrategy : public Strategy<St>
{
public:
    /*!
     * \brief simulate creates a \a Move that will be performed on the current Simulation::_state to
     * simulate a game
     * \param move the move produced by the simulation. Cast to bool to see if this is a possible
     * move.
     */
    virtual void simulate(Mv& move) noexcept = 0;

    virtual ~SimulationStrategy() noexcept = default;

protected:
private:
};

/*************************************************************************************************/
/*!
 * @brief TODO
 */
template<class St,
         class Mv,
         typename = std::enable_if_t<std::is_base_of_v<State, St>>,
         typename = std::enable_if_t<std::is_base_of_v<Move, Mv>>>
class BackPropagateStrategy : public Strategy<St>
{
public:
    virtual void backPropagate() noexcept = 0;

    virtual ~BackPropagateStrategy() noexcept = default;

protected:
private:
};

} // namespace mcts

#endif // MCTS_HPP
