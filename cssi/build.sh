#!/bin/sh

echo "rm ./obj/* ./bin/*"
rm ./obj/*
rm ./bin/*

echo ""
echo "compile naive"
make --file=makefile_NAIVE
mv bin/main bin/naive

echo ""
echo "compile lambda"
make --file=makefile_LAMBDA
mv bin/main bin/lambda

echo ""
echo "compile dfs"
make --file=makefile_DFS
mv bin/main bin/dfs
