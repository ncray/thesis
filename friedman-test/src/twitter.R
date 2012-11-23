library(twitteR)
library(plyr)
## Get new data from Twitter
## obama <- userTimeline("BarackObama", n = 1500)
## romney <- userTimeline("MittRomney", n = 1500)
## palin <- userTimeline("SarahPalinUSA", n = 1500)
## save(obama, "obama")
## save(romney, "romney")
## save(palin, "palin")

processTweet <- function(tweet){
  text <- tweet$getText()
  strippedURLs <- gsub('\\S*\\.\\S*\\/\\S*', "", text, perl = TRUE)
  newLinesAreSpaces <- gsub('\\n', " ", strippedURLs, perl = TRUE)
  extractedAlphabet <- gsub("[^ a-zA-Z]", "", newLinesAreSpaces, perl = TRUE)
  strippedWhitespaces <- gsub('^\\s+|\\s+$', "", extractedAlphabet, perl = TRUE)
  collapsedSpaces <- gsub('\\s+', " ", strippedWhitespaces, perl = TRUE)
  tolower(collapsedSpaces)
}

filterTweets <- function(txt.l){
  long.ind <- which(laply(txt.l, nchar) > 10)
  longerTweets <- txt.l[long.ind]
  unique(longerTweets)
}

processAndFilter <- function(tweets) filterTweets(llply(tweets, processTweet))

checkProcessing <- function(){
  load("obama")
  ## Pre/post processing
  l_ply(obama[1:10], function(tweet){
    print(tweet)
    print(processTweet(tweet))
  })
}

getObama <- function(){
  load("obama")
  processAndFilter(obama)
}

getPalin <- function(){
  load("palin")
  processAndFilter(palin)
}
