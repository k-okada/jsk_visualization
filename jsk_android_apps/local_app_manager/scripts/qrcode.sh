#!/bin/bash
HOST=$ROS_IP
PORT=11311
wget "http://chart.apis.google.com/chart?cht=qr&chs=350x350&chl=URL:+http://$HOST:$PORT/" -O /tmp/qrcode.png -q
gnome-open /tmp/qrcode.png