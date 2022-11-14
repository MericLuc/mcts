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

// TODO
// I know everything will be there, so that's a start !

#ifndef MCTS_HPP
#define MCTS_HPP

namespace mcts {} // namespace mcts

#endif // MCTS_HPP
