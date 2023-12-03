git config --global --unset http.proxy
git config --global --unset https.proxy
echo "reseted git proxy as follows:"
echo "git_porxy_http is:"$(git config --global --get http.proxy)
echo "git_porxy_https is:"$(git config --global --get https.proxy)