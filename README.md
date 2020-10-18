# Shortest-Path Algorithms Comparison

Comparing shortest-path algorithms: Dijkstra, A*, and Landmark

Simon Fraser University CMPT307 - Fall 2019

## Description
### Output
**First part**: Randomly generate 20 query pairs of the form (q1,q2). For each query:
* Determine the shortest distance in G from vertex q1 to vertex q2 using Dijkstra’s, A* and Landmark  algorithms. 
* Output the number of vertices Dijkstra’s algorithm has visited.
* Output the number of vertices A* algorithm has visited.
* Output the number of vertices Landmark algorithm has visited.

**Second part**: Determine, on an average, what are the savings if A* and Landmark algorithms are used instead of Dijkstra’s algorithm


### How to run
In terminal:
```bash
./a.out
```

### How it works
The [input data file](https://github.com/wendyhwl/Shortest-Path-Comparison/blob/main/graph1000.txt) contains a graph with 1000 vertices. Each vertex has a (x,y) coordinate where x is the latitude and y is the longitude of the vertex. The program takes in the input data and constructs a Graph class. Each of the shortest-path algorithm then takes in the Graph and a set of randomly generated vertices as parameters, and generates the results (the number of visited vertices and total distance between two vertices).

The average number of nodes visited is calculated using the sum of total number of nodes visited from the 20 different sets of vertices, dividing by 20. The result shows that, in terms of number of nodes visited, Landmark is the most efficient, whereas Dijkstra is the least efficient algorithm.

### Format
Randomly generates 20 query pairs in the form of (q1,q2). For each query, it outputs:
* The number of vertices visited using Dijkstra, A* and Landmark algorithms
* The shortest distance generated using Dijkstra, A* and Landmark algorithms

Output format:
```bash
QUERY #: q1 -> q2
1. Dijkstra:  number of nodes visited   (distance)
2. Astar:     number of nodes visited   (distance)
3. Landmark:  number of nodes visited   (distance)
```

At the end, the program outputs the average vertices visited for each of the above algorithm.
Output is in the form of:
```bash
Average
1. Dijkstra   average number of nodes visited
2. Astar      average number of nodes visited
3. Landmark   average number of nodes visited
```

## Language

C++

## Acknowledgement

Mini project created by CMPT310 Professor Binay Bhattacharya.

[SFU Academic Integrity](http://www.sfu.ca/students/academicintegrity.html)
