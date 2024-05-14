#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subtitle:: *)
(*This is an example of running a C code from Mathematica *)



(* ::Input:: *)
ParseFaceData[data_] := Module[{func, pos}, func = data[[All, 7]];
  pos = Table[Take[data[[i]], 6], {i, Length[data]}];
  pos = Partition[#, 2] & /@ pos;
  Table[{pos[[i]], func[[i]]}, {i, Length[pos]}]]


ParseVertexData[data_] := Module[{points, r1, r2},
  r1 = {1, 0};
  r2 = {1/2, Sqrt[3]/2};
  Table[
  points = r1*data[[i, 1]]+r2*data[[i, 2]]; 
  Flatten[{points, data[[i, 3]]}]
  , {i, Length[data]}]]

takeRatio[a_, b_]:=Module[{data = {}, c = a},
Do[c[[i, 3]] = a[[i, 3]]/b[[i, 3]], {i, Length[a]}];
c
]
Data = Import["Area.dat", "table"];
res = ParseFaceData[Data];
minValue = Min[res[[All, 2]]];
maxValue = Max[res[[All, 2]]];
scaledColorFunction[value_] := 
  ColorData["Rainbow"][(value - minValue)/(maxValue - minValue)];
heatMap = 
  Graphics[{scaledColorFunction[#[[2]]], Polygon[#[[1]]]} & /@ res, 
   Axes -> True];
legend = BarLegend[{"Rainbow", {minValue, maxValue}}];
fig = Legended[heatMap, legend]
Export["area.png", fig]

Data = Import["Perimeter.dat", "table"];
res = ParseFaceData[Data];
minValue = Min[res[[All, 2]]];
maxValue = Max[res[[All, 2]]];
scaledColorFunction[value_] := 
  ColorData["Rainbow"][(value - minValue)/(maxValue - minValue)];
heatMap = 
  Graphics[{scaledColorFunction[#[[2]]], Polygon[#[[1]]]} & /@ res, 
   Axes -> True];
legend = BarLegend[{"Rainbow", {minValue, maxValue}}];
fig = Legended[heatMap, legend]
Export["perimeter.png", fig]

Data = Import["Circumradius.dat", "table"];
res = ParseFaceData[Data];
minValue = Min[res[[All, 2]]];
maxValue = Max[res[[All, 2]]];
scaledColorFunction[value_] := 
  ColorData["Rainbow"][(value - minValue)/(maxValue - minValue)];
heatMap = 
  Graphics[{scaledColorFunction[#[[2]]], Polygon[#[[1]]]} & /@ res, 
   Axes -> True];
legend = BarLegend[{"Rainbow", {minValue, maxValue}}];
fig = Legended[heatMap, legend];
Export["circumradius.png", fig]

df = Import["deficitAngle.dat", "table"];
dual = Import["dualArr.dat", "table"];
ct= Import["counterterm.dat", "table"];
ratio = takeRatio[df, dual]; 
fig=ListDensityPlot[ParseVertexData[df], ColorFunction -> "Rainbow", PlotLegends -> True,
  FrameLabel -> {"x", "y"}];
Export["deficit.png", fig]

fig=ListDensityPlot[ParseVertexData[ratio], ColorFunction -> "Rainbow", PlotLegends -> True,
  FrameLabel -> {"x", "y"}];
Export["ratio.png", fig]

fig=ListDensityPlot[ParseVertexData[dual], ColorFunction -> "Rainbow", PlotLegends -> True,
  FrameLabel -> {"x", "y"}];
Export["dual.png", fig]


fig=ListDensityPlot[ParseVertexData[ct], ColorFunction -> "Rainbow", PlotLegends -> True,
  FrameLabel -> {"x", "y"}, PlotRange -> All];
Export["counterterm.png", fig]


