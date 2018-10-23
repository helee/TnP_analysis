#!/bin/bash

$ROOTSYS/bin/root -b -l <<EOF
{
    .L /Users/helee/temp/analysis/treeAnalysis.C
    treeAnalysis m("TnPTree_data_63.root");
    m.Loop();
}
EOF
