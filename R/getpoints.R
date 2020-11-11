#'@name getpoints
#'@title Extract the coordinates from Kaplan-Meier(K-M) curves by mouse-clicks
#'@description The \code{getpoints()} function extracts the coordinates from K-M curves by mouse-clicks. The K-M curves should be in the format of bitmap images(in JPEG,PNG,BMP,JPG or TIFF),
#'             and the use of .png file is highly recommended, since it can greatly shorten the processing time in R. \cr\cr
#'             In addition to the image itself, the input of the \code{getpoints()} function includes two x-coordinates (x1 and x2) and two y-coordinates to decide the location and scale of the curve. Once the
#'             image is read into R and displayed in the plots window, firstly the user need to click on the four points on the x-axis and y-axis according to the input, and in the order of (x1,x2,y1,and y2); secondly, the user need to collect the points coordinates by
#'             mouse-clicks on the curve. To get desirable estimation, we suggest collecting 80-100 points on each curve, and including the points where the survival probability drops. The output of this function is a two-column dataset of coordinates extracted from the K-M curve.
#'
#'
#'@usage getpoints(f,x1,x2,y1,y2)
#'@param f the bitmap image(in JPEG,PNG,BMP,JPG or TIFF formate) of the K-M curves. The input can be either the pathway to the image file, or the bitmap digital image itself.
#'@param x1 two points needed to decide the postion and scale of the x-axis. Here x1 is the actual x-coordinate of the right point on x-axis
#'@param x2 two points needed to decide the postion and scale of the x-axis. Here x2 is the actual x-coordinate of the left point on x-axis
#'@param y1 two points needed to decide the postion and scale of the y-axis. Here y1 is the actual y-coordinate of the lower point on y-axis
#'@param y2 two points needed to decide the postion and scale of the y-axis. Here y2 is the actual y-coordinate of the upper point on y-axis

#'@return \code{getpoints()} returns a two-column dataset of coordinates extracted from a K-M curve.\cr \cr
#' @importFrom graphics locator
#' @importFrom stats lm
#' @importFrom readbitmap read.bitmap
#' @importFrom graphics par plot.new rasterImage
#' @export
#' @examples
#' str(imgexp)
#'\donttest{
#'
#' ## Extract the coordinates from Kaplan-Meier(K-M) curves by mouse-clicks.
#' ## The K-M curve should be in the format of bitmap images. The input f should be either
#' ## the pathway to the image file, or the bitmap digital image itself.
#' ## Example: extract coordinates from the sample bitmap digital image (imgexp)
#' plot.new()
#' rasterImage(imgexp, 0, 0, 1, 1)
#' ## User need to use mouse-clicks to decide the positions of coordinates,
#' ## and the points want to extract.
#' df <- getpoints(imgexp,0,60,0,100)
#' head(df)
#' ## the extracted dataset df can be used to estimate IPD by other functions in the package
#' trisk <- Radiationdata$trisk
#' nrisk.radio <- Radiationdata$nrisk.radio
#' pre_radio <- preprocess(dat=df, trisk=trisk,
#'              nrisk=nrisk.radio,totalpts=NULL,maxy=100)
#' est_radio <- getIPD(prep=pre_radio,armID=0,tot.events=NULL)
#' }
#'
#'@references Poisot T. The digitize package: extracting numerical data from scatterplots. The R Journal. 2011 Jun 1;3(1):25-6.


getpoints <- function(f,x1,x2,y1,y2){

  ## if bitmap
  if (typeof(f)=="character")
      { lfname <- tolower(f)
       if ((strsplit(lfname,".jpeg")[[1]]==lfname) && (strsplit(lfname,".tiff")[[1]]==lfname) &&
            (strsplit(lfname,".bmp")[[1]]==lfname) && (strsplit(lfname,".png")[[1]]==lfname) &&
           (strsplit(lfname,".jpg")[[1]]==lfname))
         {stop ("This function can only process bitmap images in JPEG, PNG,BMP, or TIFF formate.")}
       img <- readbitmap::read.bitmap(f)
      } else if (typeof(f)=="double")
      {
        img <- f
      } else {
        stop ("Please double check the format of the image file.")
      }
  ## function to read the bitmap and points on x-axis and y-axis
  axispoints <- function(img){
    op <- par(mar = c(0, 0, 0, 0))
    on.exit(par(op))
    plot.new()
    rasterImage(img, 0, 0, 1, 1)
    message("You need to define the points on the x and y axis according to your input x1,x2,y1,y2. \n")
    message("Click in the order of left x-axis point (x1), right x-axis point(x2),
      lower y-axis point(y1), and upper y-axis point(y2). \n")
    x1 <- as.data.frame(locator(n = 1,type = 'p',pch = 4,col = 'blue',lwd = 2))
    x2 <- as.data.frame(locator(n = 1,type = 'p',pch = 4,col = 'blue',lwd = 2))
    y1 <- as.data.frame(locator(n = 1,type = 'p',pch = 3,col = 'red',lwd = 2))
    y2 <- as.data.frame(locator(n = 1,type = 'p',pch = 3,col = 'red',lwd = 2))
    ap <- rbind(x1,x2,y1,y2)
    return(ap)
  }

  ## function to calibrate the points to the appropriate coordinates
  calibrate <-  function(ap,data,x1,x2,y1,y2){
    x  <- ap$x[c(1,2)]
    y  <- ap$y[c(3,4)]
    cx <- lm(formula = c(x1,x2) ~ c(x))$coeff
    cy <- lm(formula = c(y1,y2) ~ c(y))$coeff
    data$x <- data$x*cx[2]+cx[1]
    data$y <- data$y*cy[2]+cy[1]
    return(as.data.frame(data))
   }

  ## take the points
  ap <- axispoints(img)
  message("Mouse click on the K-M curve to take the points of interest. The points will only be labled when you finish all the clicks.")
  takepoints <- locator(n=512,type='p',pch=1,col='orange4',lwd=1.2,cex=1.2)
  df <- calibrate(ap,takepoints,x1,x2,y1,y2)
  par()
  return(df)
}


