L=64;

AData = Import["Evan_data/CurvedTriangleupdate_64_final.dat"];
AData[[All, {3,4,5}]] = L*AData[[All, {3,4,5}]] /Norm/@(AData[[All, {3,4,5}]])


BData = Import["Evan_data/CurvedTriangleupdate_64_initial.dat"];
BData[[All, {3,4,5}]] = L*BData[[All, {3,4,5}]] /Norm/@(BData[[All, {3,4,5}]])


Diff = AData - BData;
(*PlotData = 
  Table[{AData[[All, {3, 4, 5}]][[i]], 
    Diff[[All, {3, 4, 5}]][[i]]}, {i, Length[AData]}];
makeArrow[{pos_, vec_}] := Arrow[{pos, pos + vec}];
gr = Graphics3D[{{Red,(*Vector color*)
    Arrowheads[
     0.02],(*Size of the arrowheads*)(makeArrow /@ 
      PlotData) (*Applying the function to each data point*)}, {Blue, 
    PointSize[Medium], Point[PlotData[[All, 1]]]}}, Boxed -> True, 
  Axes -> True, ViewPoint -> {0, 0, 1}];
  Export["displacement.pdf", gr];*)


returnFullPos[L_] := Module[{points, r1, r2}, r1 = {1, 0};
  r2 = {1/2, Sqrt[3]/2};
  Table[Table[r2*j + r1*i, {i, 0 , L - j}], {j, 0, L}]]
func = Norm /@ Diff[[All, {3, 4, 5}]];
pos = Flatten[returnFullPos[L], 1];
StackedData = Table[Flatten[{pos[[i]], func[[i]]}], {i, Length[func]}];
fig = ListDensityPlot[StackedData, ColorFunction -> "Rainbow", 
  PlotLegends -> True, FrameLabel -> {"x", "y"}, PlotRange -> All]
Export["magnitude.pdf", fig]


