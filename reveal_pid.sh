#!/bin/bash

pid=$$
echo $pid

sleep 0.5

# Carefull with this, if /volatile doesn't exist?!
#dd if=/dev/urandom of=/volatile/cjb/randdata/$pid\_rand.img bs=1M count=10
