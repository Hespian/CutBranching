## Building

To build the code run

```shell
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
from the top-level directory of this repository. Nested dissection branching uses a slightly modified version of [InertialFlowCutter](https://github.com/kit-algo/InertialFlowCutter) with [KaHIP](https://github.com/KaHIP/KaHIP) backend. Those submodules require TBB, OpenMP and MPI. If you do not want to use nested dissection branching, you can disable it by setting the option `USE_IFC` to  `OFF` in `CMakeLists.txt`.

## Computing a MIS

Input graphs have to be undirected and must not contain duplicate edges or self-loops. The code expects the following format: An input file contains the adjacency list of the graph. Nodes are numbered from `0` to `n-1`. The format for the first line is: `p NUM_NODES NUM_EDGES`. All following lines are of the form `u v` and represent an edge of the graph. Note that for an edge `{u,v}` of the graph the input file has to contain both `u v` and `v u`. Lines beginning with c are ignored 

Run

```shell
./MIS <graph_dir_path> <output_path> [<strat>]
```
to compute a MIS for all graphs in `graph_dir_path`. The result is written to the `output_path`directory. The other parameters are optional and control which strategy is used and how it is parameterized. Details on the different strategies can be found in our paper. By default, max. deg. branching is used. Viable option for the strategies are: 

- `0` for branching on a random vertex
- `1` for branching on a min. deg. vertex
- `2` for branching on a max. deg. vertex
- `3` for branching using articulation points
- `4` for branching using edge cuts
- `5` for branching using nested dissection
- `6` for twin-reduction-based branching
- `7` for funnel-reduction-based branching
- `8` for unconfined-reduction-based branching
- `9` for combined branching

## Paper

If you use this code in your publication, please cite our paper:

```
@InProceedings{hespe_et_al:LIPIcs.SEA.2021.17,
  author =	{Hespe, Demian and Lamm, Sebastian and Schorr, Christian},
  title =	{{Targeted Branching for the Maximum Independent Set Problem}},
  booktitle =	{19th International Symposium on Experimental Algorithms (SEA 2021)},
  pages =	{17:1--17:21},
  series =	{Leibniz International Proceedings in Informatics (LIPIcs)},
  ISBN =	{978-3-95977-185-6},
  ISSN =	{1868-8969},
  year =	{2021},
  volume =	{190},
  editor =	{Coudert, David and Natale, Emanuele},
  publisher =	{Schloss Dagstuhl -- Leibniz-Zentrum f{\"u}r Informatik},
  address =	{Dagstuhl, Germany},
  URL =		{https://drops.dagstuhl.de/opus/volltexte/2021/13789},
  URN =		{urn:nbn:de:0030-drops-137893},
  doi =		{10.4230/LIPIcs.SEA.2021.17},
  annote =	{Keywords: Graphs, Combinatorial Optimization, Independent Set, Vertex Cover, Clique, Branch-and-Reduce, Branch-and-Bound, Data Reduction}
}
```
