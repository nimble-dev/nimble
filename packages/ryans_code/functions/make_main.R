make_main=function(allModels){
  html="<!DOCTYPE html PUBLIC>
  <html>
  <head>
  <body>
  <h1>Examples</h1>"
  for (i in 1:length(allModels)){
    html=paste(html,"<a href=\"",allModels[i],".html\">",allModels[i],"</a>\n <br>\n",sep="")
  }
  html=paste(html,"</body>
  </head>
  </html>",sep="")
  cat(html,file='main.html')
}