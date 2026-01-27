Every time fedra gets updated
From the lxplus master branch git-svn
-- git svn rebase
-- git push
From the fedradoxygen branch

-- git pull
-- doxygen doxygen/fedoxy.conf
-- git commit -m 'update commit'
-- git rebase origin/master
-- git push