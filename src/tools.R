####### Various tools for simulations and plotting



makeTransparent<-function(someColor, alpha=70)
{ # from https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
