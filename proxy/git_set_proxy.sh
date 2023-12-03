hostip=$(cat /etc/resolv.conf | grep nameserver | awk '{ print $2 }')
wslip=$(hostname -I | awk '{print $1}')
echo "host ip is:"$hostip
echo "wsl ip is:"$wslip
port=21882
echo "port is:"$port
export https_proxy=http://${hostip}:${port} 
export http_proxy=http://${hostip}:${port} 
export all_proxy=socks5://${hostip}:${port}
echo "all proxy is config successfully!"
echo "initial git proxy as follows:"
echo "git_porxy_http is:"$(git config --global --get http.proxy)
echo "git_porxy_https is:"$(git config --global --get https.proxy)
git config --global http.proxy "http://${hostip}:${port}"
git config --global https.proxy "http://${hostip}:${port}"
echo "new git proxy as follows:"
echo "git_porxy_http is:"$(git config --global --get http.proxy)
echo "git_porxy_https is:"$(git config --global --get https.proxy)