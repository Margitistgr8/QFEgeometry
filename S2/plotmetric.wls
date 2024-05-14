AData = Import[
   "/Users/jinyunlin/QFE-Research/affineFlow-2/JintestCode/testing/\
ConstantLattice/Jin_data/CurvedTriangleupdate_64_final.dat", "Table"];
BData = Import[
   "/Users/jinyunlin/QFE-Research/affineFlow-2/JintestCode/testing/\
ConstantLattice/Evan_data/CurvedTriangleupdate_64_final.dat", "Table"];
L=64;
AData *= L/Norm@AData[[1,3;;5]];
BData *= L/Norm@BData[[1,3;;5]];




extractRowLinks[arr_] := Module[{L = Length[arr], start, end},
  start = arr[[1 ;; -2]];
  end = arr[[2 ;;]];
  Norm /@ (end - start)]
ParseFaceData[data_] := Module[{r1, r2}, 
  r1 = {1, 0};
  r2 = {1/2, Sqrt[3]/2};
  Table[Flatten[{data[[i, 1]]*r1 + r2*data[[i, 2]], 
     data[[i, 3]]}], {i, Length[data]}]
  ]
  returnstruct[data_, L_] := Module[{counter = L, ini = 1, arr},
  arr = {} ;
  Do[AppendTo[arr, data[[ini ;; ini + counter]]];
   ini += counter + 1; 
   counter -= 1;
   , {i, 1, L}]; 
  arr
  ]
returnPos[L_] := Module[{points, r1, r2},
  r1 = {1, 0};
  r2 = {1/2, Sqrt[3]/2};
  Table[Table[r2*j + r1*i, {i, 0, L - j - 1}], {j, 0, L - 1}]]


res = returnstruct[AData, L];
edgesA = extractRowLinks /@ res;
res = returnstruct[BData, L];
edgesB = extractRowLinks /@ res;

pos = returnPos[L];
edges = edgesB-edgesA; 
Data = Table[
   Flatten[{Flatten[pos, 1][[i]], Flatten[edges, 1][[i]]}], {i, 
    Length[Flatten[pos, 1]]}]; 

fig =ListDensityPlot[Data, ColorFunction -> "Rainbow", PlotLegends -> True,
  FrameLabel -> {"x", "y"}]; 

Export["Metric_Diff.pdf", fig]