#include <assert.h>
#include <stdio.h>
#include <sys/wait.h>
#include <unistd.h>
#include <iostream>
#include <vector>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define ABS(a) ((a) < 0 ? (-(a)) : (a))

typedef int number;

struct ConfigurationView {
  number size;    // size of square box
  number** grid;  // built to be two units larger than 'size'. Should be indexed
                  // from 1 to 'size' (index 0 and size - 1 will contain -1 to
                  // indicate the box)
};

ConfigurationView* createConfigurationView(number size) {
  auto conf = new ConfigurationView;
  conf->size = size;
  conf->grid = new number*[size + 2];
  for (number i = 0; i < size + 2; i++) {
    conf->grid[i] = new number[size + 2];

    // fill the box edges with -1
    if (i == 0 || i == size + 1) {
      for (number col = 0; col < size + 2; col++) {
        conf->grid[i][col] = -1;
      }
    } else {
      conf->grid[i][0] = conf->grid[i][size + 1] = -1;
    }
  }

  return conf;
}

ConfigurationView* copyConfigurationView(ConfigurationView* conf,
                                         ConfigurationView* target) {
  // delete old grid, if any
  if (target->grid) {
    for (number i = 0; i < target->size + 2; i++) {
      delete[] target->grid[i];
    }
    delete[] target->grid;
  }

  // copy grid
  target->size = conf->size;
  target->grid = new number*[target->size + 2];
  for (number i = 0; i < target->size + 2; i++) {
    target->grid[i] = new number[target->size + 2];

    for (number col = 0; col < target->size + 2; col++) {
      target->grid[i][col] = conf->grid[i][col];
    }
  }
}

void deleteConfigurationView(ConfigurationView* conf) {
  // delete grid, if any
  if (conf->grid) {
    for (number i = 0; i < conf->size + 2; i++) {
      delete[] conf->grid[i];
    }
    delete[] conf->grid;
  }
  delete conf;
}

void printConfigurationView(ConfigurationView* conf) {
  using namespace std;
  cout << "ConfigurationView with size " << conf->size << endl;
  for (number y = conf->size + 1; y >= 0; y--) {
    for (number x = 0; x < conf->size + 2; x++) {
      printf("%3d", conf->grid[x][y]);
    }
    cout << endl;
  }

  cout << "End grid" << endl;
}

struct Square {
  number size;
  number area;
  bool packed;
};

struct Neighbours {
  number left, top, right, bottom;
};

struct Metrics {
  float cavingDegree;
  number cornerDegree;
  number edgeDegree;
  number backtrackingValue = 0;
};

struct Packing {
  number x;
  number y;
  number size;
  number area;
  Square* square;
  Neighbours neighbours;
  Metrics metrics;
};

template <class T>
struct List {
  T current;

  // these two should not be relied upon for detecting end-of-list. Instead, use
  // the 'length' and the 'index' attributes
  List* next = nullptr;
  List* before = nullptr;

  // valid for all nodes
  List* head = this;
  // only valid for the head node
  number length = 0;
  // valid for all nodes
  number index = 0;
  // valid only for the head
  // TODO: Currently, not valid for any node. DO NOT USE
  List* last = this;
  inline List* getNext(bool forceLengthIncrease = false) {
    if (next && forceLengthIncrease) {
      head->length++;
      head->last = next;
    }
    if (!next) {
      next = new List;
      next->next = nullptr;
      next->before = this;
      head->length++;
      next->index = index + 1;
      next->head = head;
      head->last = next;
    }
    return next;
  }
  /**
   * Makes it look like there is no node in the list (except the head node,
   * which is always there)
   */
  inline void softReset() {
    head->length = 0;
    head->last = head;
  }
  /**
   * Makes it look like the last node is not there anymore
   */
  inline void softDecrement() {
    head->length--;
    head->last = head->last->before;
  }
  /**
   * Makes it look like 'node' is the last node in the list.
   *
   * If 'node' was not originally a node in the list, things will go terribly
   * wrong.
   */
  inline void softResetTo(List* node) {
    head->length = node->index;
    head->last = node;
  }
};

#define FOREACH(init, list, action)                       \
  {                                                       \
    init = (list);                                        \
    while (true) {                                        \
      init = init->next;                                  \
      if (!init || (init->index > (list)->length)) break; \
      {action};                                           \
    };                                                    \
  }

template <class T>
inline List<T>* beforeLast(List<T>& list) {}

struct Configuration {
  number boxSize;
  // should use a lead node
  List<Packing>* packings;
};

// manhattan distance for rectangles
inline number distance(const Packing& p1, const Packing& p2) {
  return MAX(p1.x < p2.x ? (p2.x - (p1.x + p1.size))
                         : (p1.x - (p2.x + p2.size)),
             0) +
         MAX(p1.y < p2.y ? (p2.y - (p1.y + p1.size))
                         : (p1.y - (p2.y + p2.size)),
             0);
  // return MAX(0, (ABS(p1.x - p2.x) << 1) - (p1.size + p2.size)) +
  //        MAX(0, (ABS(p1.y - p2.y) << 1) - (p1.size + p2.size));
}

// zero if on box edge
inline number distanceToBox(const Packing& p, number boxSize) {
  return MIN(MIN(p.x, boxSize - p.x - p.size),
             MIN(p.y, boxSize - p.y - p.size));
}

// considering only non-zero distances to box walls
inline number minDistToWalls(const Packing& p, number boxSize) {
  number dist = boxSize;
  if (p.x != 0) dist = MIN(dist, p.x);
  if (p.y != 0) dist = MIN(dist, p.y);
  if (p.x + p.size != boxSize) dist = MIN(dist, boxSize - (p.x + p.size));
  if (p.y + p.size != boxSize) dist = MIN(dist, boxSize - (p.y + p.size));
  return dist;
}

inline bool onBoxCorner(const Packing& p, number boxSize) {
  return MIN(p.x, boxSize - p.x - p.size) == 0 &&
         MIN(p.y, boxSize - p.y - p.size) == 0;
}

inline bool overlap(const Packing& p1, const Packing& p2) {
  return (p1.x < p2.x ? p2.x < p1.x + p1.size : p1.x < p2.x + p2.size) &&
         (p1.y < p2.y ? p2.y < p1.y + p1.size : p1.y < p2.y + p2.size);
};

// true iff 'p1' overlaps any packing in the list
inline bool overlap(Packing& p1, List<Packing>* list) {
  List<Packing>* pack;
  FOREACH(pack, list, if (overlap(p1, pack->current)) { return true; });
  return false;
};

inline bool insideBox(const Packing& p, number boxSide) {
  return p.x >= 0 && p.x + p.size <= boxSide && p.y >= 0 &&
         p.y + p.size <= boxSide;
};

// true iff 'p1' and 'p2' share an edge
inline bool neighbour(const Packing& p1, const Packing& p2) {
  return distance(p1, p2) == 0 &&
         !((p1.x + p1.size == p2.x || p2.x + p2.size == p1.x) &&
           (p1.y + p1.size == p2.y || p2.y + p2.size == p1.y));
}

// how many left,right,bottom,top neighbours 'coa' has in 'config'
inline void countNeighbours(Configuration& config, const Packing& coa,
                            Neighbours& target) {
  target.bottom = 0;
  target.left = 0;
  target.right = 0;
  target.top = 0;
  List<Packing>* p;
  FOREACH(p, config.packings, {
    if (!neighbour(coa, p->current)) continue;
    if (coa.x + coa.size == p->current.x)
      target.right++;
    else if (p->current.x + p->current.size == coa.x)
      target.left++;
    else if (coa.y + coa.size == p->current.y)
      target.top++;
    else if (p->current.y + p->current.size == coa.y)
      target.bottom++;
  });
  // check box edges
  if (coa.x == 0) target.left++;
  if (coa.x + coa.size == config.boxSize) target.right++;
  if (coa.y == 0) target.bottom++;
  if (coa.y + coa.size == config.boxSize) target.top++;
}

inline bool onCorner(Neighbours& neighbours) {
  return (neighbours.left > 0 || neighbours.right > 0) &&
         (neighbours.bottom > 0 || neighbours.top > 0);
}

void enumerateCOAs(Configuration& config, List<Square>& squares,
                   List<Packing>& target) {
  List<Square>* sq;
  List<Packing>* pivot;
  target.length = 0;  // critical !
  auto packing = target.getNext(true);
  FOREACH(sq, &squares, {
    if (sq->current.packed) continue;
    packing->current.size = sq->current.size;
    packing->current.area = sq->current.area;
    packing->current.square = &sq->current;
    for (number x = 0; x <= config.boxSize - sq->current.size; x++) {
      for (number y : {0, config.boxSize - sq->current.size}) {
        packing->current.x = x;
        packing->current.y = y;
        if (overlap(packing->current, config.packings)) continue;
        countNeighbours(config, packing->current, packing->current.neighbours);
        if (!onCorner(packing->current.neighbours)) continue;
        // prepare next node
        packing = packing->getNext(true);
        packing->current.size = sq->current.size;
        packing->current.area = sq->current.area;
        packing->current.square = &sq->current;
      }
    }
    for (number y = 0; y <= config.boxSize - sq->current.size; y++) {
      for (number x : {0, config.boxSize - sq->current.size}) {
        packing->current.x = x;
        packing->current.y = y;
        if (overlap(packing->current, config.packings)) continue;
        countNeighbours(config, packing->current, packing->current.neighbours);
        if (!onCorner(packing->current.neighbours)) continue;
        // prepare next node
        packing = packing->getNext(true);
        packing->current.size = sq->current.size;
        packing->current.area = sq->current.area;
        packing->current.square = &sq->current;
      }
    }
    FOREACH(pivot, config.packings, {
      // pack under the pivot
      packing->current.y = pivot->current.y - sq->current.size;
      if (packing->current.y >= 0) {
        for (number x = MAX(0, pivot->current.x - sq->current.size + 1);
             x <= MIN(config.boxSize - sq->current.size,
                      pivot->current.x + pivot->current.size - 1);
             x++) {
          packing->current.y = pivot->current.y - sq->current.size;
          packing->current.x = x;
          if (overlap(packing->current, config.packings)) continue;
          countNeighbours(config, packing->current,
                          packing->current.neighbours);
          if (!onCorner(packing->current.neighbours)) continue;
          // prepare next node
          packing = packing->getNext(true);
          packing->current.size = sq->current.size;
          packing->current.area = sq->current.area;
          packing->current.square = &sq->current;
        }
      }
      // pack above the pivot
      packing->current.y = pivot->current.y + pivot->current.size;
      if (packing->current.y + packing->current.size <= config.boxSize) {
        for (number x = MAX(0, pivot->current.x - sq->current.size + 1);
             x <= MIN(config.boxSize - sq->current.size,
                      pivot->current.x + pivot->current.size - 1);
             x++) {
          packing->current.y = pivot->current.y + pivot->current.size;
          packing->current.x = x;
          if (overlap(packing->current, config.packings)) continue;
          countNeighbours(config, packing->current,
                          packing->current.neighbours);
          if (!onCorner(packing->current.neighbours)) continue;
          // prepare next node
          packing = packing->getNext(true);
          packing->current.size = sq->current.size;
          packing->current.area = sq->current.area;
          packing->current.square = &sq->current;
        }
      }
      // pack to the left of the pivot
      packing->current.x = pivot->current.x - sq->current.size;
      if (packing->current.x >= 0) {
        for (number y = MAX(0, pivot->current.y - sq->current.size + 1);
             y <= MIN(config.boxSize - sq->current.size,
                      pivot->current.y + pivot->current.size - 1);
             y++) {
          packing->current.x = pivot->current.x - sq->current.size;
          packing->current.y = y;
          packing->current.x = pivot->current.x - sq->current.size;
          packing->current.y = y;
          if (overlap(packing->current, config.packings)) continue;
          countNeighbours(config, packing->current,
                          packing->current.neighbours);
          if (!onCorner(packing->current.neighbours)) continue;
          // prepare next node
          packing = packing->getNext(true);
          packing->current.size = sq->current.size;
          packing->current.area = sq->current.area;
          packing->current.square = &sq->current;
        }
      }
      // pack to the right of the pivot
      packing->current.x = pivot->current.x + pivot->current.size;
      if (packing->current.x + packing->current.size <= config.boxSize) {
        for (number y = MAX(0, pivot->current.y - sq->current.size + 1);
             y <= MIN(config.boxSize - sq->current.size,
                      pivot->current.y + pivot->current.size - 1);
             y++) {
          packing->current.x = pivot->current.x + pivot->current.size;
          packing->current.y = y;
          if (overlap(packing->current, config.packings)) continue;
          countNeighbours(config, packing->current,
                          packing->current.neighbours);
          if (!onCorner(packing->current.neighbours)) continue;
          // prepare next node
          packing = packing->getNext(true);
          packing->current.size = sq->current.size;
          packing->current.area = sq->current.area;
          packing->current.square = &sq->current;
        }
      }
    })
  });
  target.softDecrement();
};

// alternative caving degree definition
float cavingDegree(Configuration& config, Packing& coa) {
  number neighbours = 0;
  if (onBoxCorner(coa, config.boxSize))
    neighbours = 2;
  else if (distanceToBox(coa, config.boxSize) == 0)
    neighbours = 1;
  number minDist = minDistToWalls(coa, config.boxSize);
  List<Packing>* p;
  FOREACH(p, config.packings, number dist = distance(coa, p->current);
          if (dist == 0) {
            neighbours++;
            if (neighbours == 3) return 1;
          } else minDist = MIN(minDist, dist););
  return 1 - ((float)minDist) / ((float)coa.size);
}

// assumes Packing::Neighbours is already filled
inline number edgeDegree(Packing& coa) {
  return coa.neighbours.bottom + coa.neighbours.left + coa.neighbours.top +
         coa.neighbours.right;
}

// assumes Packing::Neighbours is already filled
inline number cornerDegree(Packing& coa) {
  return (coa.neighbours.left != 0 && coa.neighbours.bottom != 0 ? 1 : 0) +
         (coa.neighbours.left != 0 && coa.neighbours.top != 0 ? 1 : 0) +
         (coa.neighbours.right != 0 && coa.neighbours.bottom != 0 ? 1 : 0) +
         (coa.neighbours.right != 0 && coa.neighbours.top != 0 ? 1 : 0);
}

void printConfiguration(Configuration* conf, ConfigurationView* view) {
  List<Packing>* pack;
  for (number x = 1; x < view->size + 1; x++) {
    for (number y = 1; y < view->size + 1; y++) {
      view->grid[x][y] = 0;
    }
  }

  FOREACH(pack, conf->packings, {
    for (number x = pack->current.x + 1;
         x <= pack->current.x + pack->current.size; x++) {
      for (number y = pack->current.y + 1;
           y <= pack->current.y + pack->current.size; y++) {
        view->grid[x][y] = pack->current.size;
      }
    };
  });
};

void printCOA(Packing coa, ConfigurationView* view) {
  for (number x = coa.x + 1; x <= coa.x + coa.size; x++) {
    for (number y = coa.y + 1; y <= coa.y + coa.size; y++) {
      view->grid[x][y] = coa.size;
    }
  }
}

/**
 * Assumes all COAs have Packing::Neighbour already filled.
 *
 * Assumes all Metrics::backtrackingValue are already filled.
 */
Packing* chooseBestCOA(Configuration* config, List<Packing>& coas) {
  if (coas.length == 0) {
    return nullptr;
  }
  List<Packing>* coa;
  coa = coas.next;
  List<Packing>* bestCoa = coa;
  bestCoa->current.metrics.cavingDegree = cavingDegree(*config, coa->current);
  bestCoa->current.metrics.cornerDegree = cornerDegree(coa->current);
  bestCoa->current.metrics.edgeDegree = edgeDegree(coa->current);
  // we iterate the list from lower index to greater index,
  // hence we don't need to select based on index when all other metrics
  // are the same
  FOREACH(coa, &coas, {
    if (coa->current.metrics.backtrackingValue >
        bestCoa->current.metrics.backtrackingValue) {
      bestCoa = coa;
      bestCoa->current.metrics.cavingDegree =
          cavingDegree(*config, coa->current);
      bestCoa->current.metrics.cornerDegree = cornerDegree(coa->current);
      bestCoa->current.metrics.edgeDegree = edgeDegree(coa->current);
    } else if (coa->current.metrics.backtrackingValue ==
               bestCoa->current.metrics.backtrackingValue) {
      coa->current.metrics.cavingDegree = cavingDegree(*config, coa->current);
      if (coa->current.metrics.cavingDegree >
          bestCoa->current.metrics.cavingDegree) {
        bestCoa = coa;
        bestCoa->current.metrics.cornerDegree = cornerDegree(coa->current);
        bestCoa->current.metrics.edgeDegree = edgeDegree(coa->current);
      } else if (coa->current.metrics.cavingDegree ==
                 bestCoa->current.metrics.cavingDegree) {
        coa->current.metrics.cornerDegree = cornerDegree(coa->current);
        if (coa->current.metrics.cornerDegree >
            bestCoa->current.metrics.cornerDegree) {
          bestCoa = coa;
          bestCoa->current.metrics.edgeDegree = edgeDegree(coa->current);
        } else if (coa->current.metrics.cornerDegree ==
                   bestCoa->current.metrics.cornerDegree) {
          coa->current.metrics.edgeDegree = edgeDegree(coa->current);
          if (coa->current.metrics.edgeDegree >
              bestCoa->current.metrics.edgeDegree) {
            bestCoa = coa;
          } else if (coa->current.metrics.edgeDegree ==
                     bestCoa->current.metrics.edgeDegree) {
            if ((coa->current.x < bestCoa->current.x) ||
                ((coa->current.x == bestCoa->current.x) &&
                 (coa->current.y < bestCoa->current.y))) {
              bestCoa = coa;
            }
          }
        }
      }
    }
  })

  return &bestCoa->current;
}

void printSolution(Configuration& conf) {
  using namespace std;
  List<Packing>* curr;
  number count = 0;
  number area = 0;
  auto view = createConfigurationView(conf.boxSize);
  printConfiguration(&conf, view);
  FOREACH(curr, conf.packings, {
    printCOA(curr->current, view);
    count++;
    area += curr->current.area;
  });
  printConfigurationView(view);
  std::cout << "quadrados: " << count << endl;
  std::cout << "area ocupada: " << area << endl;
  deleteConfigurationView(view);
}

/**
 * Auxiliary data structures to run backtracking algorith.
 *
 * Using these avoids reallocating them at every bracktracking level.
 */
struct BacktrackHelpers {
  // should not be used by clients of backtracking algorithms. Only by the
  // algorithms themselves
  List<Packing> branchCoas;
  // should not be used by clients of backtracking algorithms. Only by the
  // algorithms themselves
  List<Packing> coas;
  // client should set
  List<Packing>* lastConfPacking;
  // Assigned by the algorithm
  number totalArea = 0;
  // client should set

  bool useSubprocesses = false;
};

/**
 * Modifies conf.packings with the best complete packing we could find.
 */
void onewayBacktrack(Configuration& conf, List<Square>& squares,
                     BacktrackHelpers* h) {
  using namespace std;

  List<Packing>& branchCoas = h->branchCoas;
  List<Packing>* currBranchCoa;
  List<Packing>& coas = h->coas;
  auto pack = h->lastConfPacking;
  number index = pack->index;
  auto packBackup = pack;
  number totalArea = h->totalArea;
  while (true) {
    // std::cout << "index now: " << index << endl;
    branchCoas.softReset();
    enumerateCOAs(conf, squares, branchCoas);
    FOREACH(currBranchCoa, &branchCoas, {
      pack = packBackup;
      number area = totalArea + currBranchCoa->current.area;
      pack = pack->getNext(true);
      pack->current = currBranchCoa->current;
      currBranchCoa->current.square->packed = true;
      while (true) {
        coas.softReset();
        enumerateCOAs(conf, squares, coas);
        Packing* bestCoa = chooseBestCOA(&conf, coas);
        if (!bestCoa) break;
        bestCoa->square->packed = true;
        pack = pack->getNext(true);
        pack->current = *bestCoa;
        area += pack->current.area;
      }
      FOREACH(pack, conf.packings, {
        if (pack->index > index) pack->current.square->packed = false;
      });
      currBranchCoa->current.metrics.backtrackingValue = area;
      // std::cout << "trying, found area: " << area << endl;
      // TODO: This is too slow !!!! (dont know why -.-)
      // conf.packings->softResetTo(packBackup);
      conf.packings->length = index;
    });
    Packing* bestBranchCoa = chooseBestCOA(&conf, branchCoas);
    if (!bestBranchCoa) break;
    totalArea += bestBranchCoa->area;
    packBackup = packBackup->getNext(true);
    packBackup->current = *bestBranchCoa;
    packBackup->current.square->packed = true;
    index++;
  }
  h->totalArea = totalArea;
  // std::cout << "total area found in 1-way backtrack: " << h->totalArea <<
  // endl;
}

void nwayBacktrack(Configuration& conf, List<Square>& squares,
                   std::vector<BacktrackHelpers*>& helpers, unsigned int n) {
  if (n == 1) {
    // std::cout << n << " base " << std::endl;
    onewayBacktrack(conf, squares, helpers[0]);
    return;
  }

  List<Packing>* lastPackBackup = helpers[n - 1]->lastConfPacking;
  int index = lastPackBackup->index;

  List<Packing>& branchCoas = helpers[n - 1]->branchCoas;
  List<Packing>* currBranchCoa;

  List<Packing>* pack;
  enumerateCOAs(conf, squares, branchCoas);
  std::cout << n << " branch size " << branchCoas.length << std::endl;

  if (!helpers[n - 1]->useSubprocesses) {
    FOREACH(currBranchCoa, &branchCoas, {
      conf.packings->length = index;
      helpers[n - 2]->totalArea =
          helpers[n - 1]->totalArea + currBranchCoa->current.area;
      pack = lastPackBackup->getNext(true);
      pack->current = currBranchCoa->current;
      helpers[n - 2]->lastConfPacking = pack;
      currBranchCoa->current.square->packed = true;
      if (currBranchCoa->index % 100 == 0) {
        std::cout << n << " FOREACH " << n - 1 << " backtrack" << std::endl;
      }
      nwayBacktrack(conf, squares, helpers, n - 1);
      FOREACH(pack, conf.packings, {
        if (pack->index > index) pack->current.square->packed = false;
      });
      // really ?
      currBranchCoa->current.metrics.backtrackingValue =
          helpers[n - 2]->totalArea;
      // index++;
    });
  } else {
    int limit = 2;
    int closed_until = -1;
    int* pipes = new int[branchCoas.length << 1];
    int* results = new int[branchCoas.length << 1];  // we use only even
    FOREACH(currBranchCoa, &branchCoas, {
      if (pipe(&pipes[(currBranchCoa->index - 1) << 1]) == -1) {
        std::cout << "aborted on " << currBranchCoa->index << std::endl;
        perror("");
        abort();
      }
      auto amIParent = (fork() != 0);
      if (!amIParent) {
        conf.packings->length = index;
        helpers[n - 2]->totalArea =
            helpers[n - 1]->totalArea + currBranchCoa->current.area;
        pack = lastPackBackup->getNext(true);
        pack->current = currBranchCoa->current;
        helpers[n - 2]->lastConfPacking = pack;
        currBranchCoa->current.square->packed = true;
        if (currBranchCoa->index % 100 == 0) {
          std::cout << n << " FOREACH " << n - 1 << " backtrack" << std::endl;
        }
        nwayBacktrack(conf, squares, helpers, n - 1);
        FOREACH(pack, conf.packings, {
          if (pack->index > index) pack->current.square->packed = false;
        });
        // really ?
        currBranchCoa->current.metrics.backtrackingValue =
            helpers[n - 2]->totalArea;
        // index++;
        write(pipes[((currBranchCoa->index - 1) << 1) + 1],
              &currBranchCoa->current.metrics.backtrackingValue,
              sizeof(currBranchCoa->current.metrics.backtrackingValue));
        // int test = 12;
        // write(pipes[((currBranchCoa->index - 1) << 1) + 1], &test,
        //       sizeof(test));
        exit(0);
      }

      if (currBranchCoa->index % limit == 0) {
        for (size_t i = closed_until + 1; i < (currBranchCoa->index << 1);
             i += 2) {
          read(pipes[i], &results[i], sizeof(int));
          close(pipes[i]);
          close(pipes[i + 1]);
        }
        closed_until = (currBranchCoa->index << 1) - 1;
        while (wait(NULL) > 0) {
        }
        std::cout << "nWayBacktrack(conf.length=" << conf.packings->length
                  << ", n==2) completed " << ((closed_until + 1) >> 1)
                  << " subprocesses " << std::endl;
      }
    });

    for (size_t i = closed_until + 1; i < branchCoas.length << 1; i += 2) {
      read(pipes[i], &results[i], sizeof(int));
      close(pipes[i]);
      close(pipes[i + 1]);
      while (wait(NULL) > 0) {
      }
    }

    FOREACH(currBranchCoa, &branchCoas, {
      currBranchCoa->current.metrics.backtrackingValue =
          results[(currBranchCoa->index - 1) << 1];
      // assert(currBranchCoa->current.metrics.backtrackingValue == 12);
    });

    delete[] pipes;
    delete[] results;

    std::cout << "asserted all" << std::endl;
  }

  conf.packings->length = index;
  Packing* bestBranchCoa = chooseBestCOA(&conf, branchCoas);
  if (!bestBranchCoa) return;
  helpers[n - 1]->totalArea = helpers[n - 1]->totalArea + bestBranchCoa->area;
  lastPackBackup = lastPackBackup->getNext(true);
  lastPackBackup->current = *bestBranchCoa;
  lastPackBackup->current.square->packed = true;
  helpers[n - 1]->lastConfPacking = lastPackBackup;
  std::cout << n << " END " << n - 1 << " backtrack" << std::endl;
  std::cout << "conf.packing now with length " << conf.packings->length
            << std::endl;
  nwayBacktrack(conf, squares, helpers, n);
}

/**
 * REMEMBER TO ZERO THE LENGTH
 */
int main(void) {
  using namespace std;

  List<Square> squares;
  auto sq = &squares;
  for (number size = 78; size > 0; size--) {
    sq = sq->getNext(true);
    sq->current.size = size;
    sq->current.area = size * size;
    sq->current.packed = false;
  }

  Configuration conf;
  conf.boxSize = 100;
  conf.packings = new List<Packing>;

  BacktrackHelpers helpers;
  helpers.lastConfPacking = conf.packings;

  // auto p = conf.packings->getNext(true);
  // p->current.area = 1;
  // p->current.size = 1;
  // p->current.x = 0;
  // p->current.y = 0;
  // p->current.square = &sq->current;
  // p->current.square->packed = true;

  // helpers.lastConfPacking = p;
  // helpers.totalArea = 1;

  unsigned int n = 2;
  vector<BacktrackHelpers*> helpersVector;
  for (size_t i = 0; i < n; i++) {
    BacktrackHelpers* h = new BacktrackHelpers;
    helpersVector.push_back(h);
  }

  // onewayBacktrack(conf, squares, helpers);
  helpersVector[n - 1]->totalArea = 0;
  helpersVector[n - 1]->lastConfPacking = conf.packings;
  helpersVector[n - 1]->useSubprocesses = true;
  nwayBacktrack(conf, squares, helpersVector, n);

  // cout << "here4" << endl;
  List<Packing>* curr;
  number count = 0;
  number area = 0;
  auto view = createConfigurationView(conf.boxSize);
  printConfiguration(&conf, view);
  FOREACH(curr, conf.packings, printCOA(curr->current, view); count++;
          area += curr->current.area;);
  printConfigurationView(view);
  std::cout << "quadrados: " << count << endl;
  std::cout << "area ocupada: " << area << endl;
};