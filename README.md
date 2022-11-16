# MCTS - Monte Carlo Tree Search

**:star2: MCTS - C++ framework for Monte-Carlo Tree Search :star2:**

Yet another implementation of [Monte-Carlo Tree Search algorithm](https://en.wikipedia.org/wiki/Monte_Carlo_tree_search).

## Requirements

- A c++17 compliant compiler.
- CMake > 3.10 (if you plan to use it as a library)

## How to use 

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

## TODO 

- [ ] Actually write the code...
- [ ] Provide explanations (and resources) on how Monte-Carlo Tree Search works.
- [ ] Explain how to use the mcts library (prerequisites).
- [ ] Provide an example 
- [ ] Link to my project(s) that use it. 