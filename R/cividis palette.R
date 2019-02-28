#cividis color pallete http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0199239#sec011
civ1<-read.csv('C:/Users/dmw63/Desktop/My documents h/color palettes/cividis rgb.csv',header = FALSE)
civ1.pal<-rgb(civ1[,1], civ1[,2], civ1[,3])
plot(1:256, col=civ1.pal)