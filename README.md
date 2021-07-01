# QC21su_Mario_4_QuantumGraphColoring
### The Graph Coloring Problem
Given a graph composed of nodes and edges, color the nodes using C colors such that no two nodes sharing an edge have the same color.

### Some Applications of Graph Coloring
    Scheduling - For example, we need to schedule exams in a way that no two exams with the same students occur at the same time. 
    
    Map Coloring - Coloring the regions of a map such that no two regions sharing a border have the same color
    
    Frequency Assignment - No two towers within range of each other use the same frequency
    
    Sudoku - No two cells in the same row, column, or subsquare may have the same color. 

### My QUBO equation
<img src="https://latex.codecogs.com/gif.latex?\text{min}\left(\left(\sum_{i&space;\in&space;c}&space;\left(\sum_{(j,&space;k)&space;\in&space;E}&space;(-s_jc_i-s_kc_i&plus;2s_jc_is_kc_i&plus;1)\right&space;)&space;\right&space;)&space;&plus;&space;\gamma&space;\left(&space;\sum_{i&space;\in&space;n}&space;\left(-\sum_{j&space;\in&space;c}&space;s_ic_j&space;&plus;&space;2\sum_{j&space;\in&space;c}\sum_{k&space;>&space;j&space;\in&space;c}&space;s_ic_js_ic_k&space;&plus;&space;1&space;\right&space;)&space;\right&space;)\right&space;)" title="\text{min}\left(\left(\sum_{i \in c} \left(\sum_{(j, k) \in E} (-s_jc_i-s_kc_i+2s_jc_is_kc_i+1)\right ) \right ) + \gamma \left( \sum_{i \in n} \left(-\sum_{j \in c} s_ic_j + 2\sum_{j \in c}\sum_{k > j \in c} s_ic_js_ic_k + 1 \right ) \right )\right )" />

s represents the nodes and c represents the color. Keep in mind that any pair of s and c are not separate variables, they form part of one vairable, but I set the equation in this form to make a disinction between the node and the color.

### How to run the program
D-Wave Ocean sdk and matplotlib are required. 

    1. Modify the nodes.txt and edges.txt files following the same format shown.
    2. To modify the number of colors used, go to colorGraph.py file and modify the colors variable to the value you want 
                    

