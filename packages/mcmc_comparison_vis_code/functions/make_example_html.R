make_example_html=function(allModels,time,height){
  for(i in 1:length(allModels)){
    html=paste("<!DOCTYPE html PUBLIC>
    <html>
    <head>
    <link rel='stylesheet' type='text/css' href='style.css'/>
    </head>
    <body>
    <h1>",allModels,"</h1>
    <a href='#table'>Time Statistics</a><br>
    <a href='#eff'>Effective Size</a><br>
    <a href='#post'>Posterior Statistics</a><br>
    <h2 id='table'>Time Statistics</h2>",
    print(xtable(time,digits=c(0,0,4)),
                type='html'),
    "<br>
    <img id='eff' src=\"",allModels,"_effectiveSize.jpg\" height=\"",height*100,"\" width=\"1200\"></img>
    <br>
    <img id='post' src=\"",allModels,"_rawData.jpg\" height=\"",height*100,"\" width=\"1200\"></img>
    </body>
    </html>",sep="")
    cat(html,file=paste(allModels,".html",sep=""))
  }
}