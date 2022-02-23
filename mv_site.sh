#!/bin/bash
cp ./report/report_main.html ./docs/index.html
rm -rf ./docs/report_main_files
cp -R ./report/report_main_files ./docs
