COMMENT=$1

git add .
git commit -c "$COMMENT"
git push origin main
