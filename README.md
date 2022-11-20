# MCTS - Monte Carlo Tree Search

**:star2: MCTS - C++ framework for Monte-Carlo Tree Search :star2:**

Yet another implementation of [Monte-Carlo Tree Search algorithm](https://en.wikipedia.org/wiki/Monte_Carlo_tree_search).

## Requirements

- A c++17 compliant compiler.
- CMake > 3.10 (if you plan to use it as a library)

## How to build/install

**mcts** is packaged as a library, and you can use it as such :

```
add_subdirectory(mcts EXCLUDE_FROM_ALL)
target_link_libraries(${PROJECT_NAME} PRIVATE mcts)
```

But you might find it overkill. 

Since **mcts** is a _header-only_ library, you can just add its header **mcts.hpp** to your project and include it where needed... 

### Build

```
[~/builds/mcts] cmake -S ~${YOUR_ECV_PATH} -DCMAKE_INSTALL_PREFIX=${YOUR_INSTALL_DIR}
-- Configuring done
-- Generating done
-- Build files have been written to: ${YOUR_INSTALL_DIR}/mcts
```

### Install

```
[~/builds/mcts] make install
[100%] Built target mcts
Install the project...
-- Install configuration: ""
-- Installing: ${YOUR_INSTALL_DIR}/mcts/lib/mcts.a
-- Installing: ${YOUR_INSTALL_DIR}/mcts/include/mcts.hpp
```

## How to use

In order to use **mcts**, you are required to provide implementations of the following interfaces :
  - **State** : State of the game.
  - **Move** : Move of the game.
  - **TerminalCriteria** : Criteria that allows to know when a State is final.
  - **TerminalEval** : Evaluation function for a terminal State.
  - **ExpansionStrategy** : Controls how the search tree is expanded.
  - **SimulationStrategy** : Allows to self-play a game until it is final.
  - **BackPropagationStrategy** : Controls the back propagation of the evaluation on the search tree.

You might want to look at [**mcts.hpp**](./mcts.hpp) for more informations.

## Applications

### Tic Tac Toe

- Sources : https://github.com/MericLuc/mcts-tic-tac-toe
- Online demo : https://mericluc.github.io/mcts/tic-tac-toe/app.html

## TODO 

- [x] Actually write the code...
- [ ] Remove stupid Node implementation (We do not need pointers)
- [ ] Debug the "MCTS::advance()" function

## Acknoledgements

The whole project is fully-based on [this paper](https://dke.maastrichtuniversity.nl/m.winands/documents/pMCTS.pdf). 

Another amazingly helpful article on [MCTS methods](http://www.incompleteideas.net/609%20dropbox/other%20readings%20and%20resources/MCTS-survey.pdf)

You might want to check **better implementations** on which it is largely inspired :
- https://github.com/Konijnendijk/cpp-mcts
- https://github.com/steve3003/mcts-cpp