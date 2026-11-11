#### SpaMTP MALDI to Visium Spot Merging Functions  ####################################################################################################################################################################


#' This function computes the coordinates of a square's four corners based on a
#' given center point and width.
#'
#' @param centerX Numeric value defining the x-coordinate of the square's
#'   center.
#' @param centerY Numeric value defining the y-coordinate of the square's
#'   center.
#' @param width Numeric value indicating the width of the square.
#' @param name Character string specifying the label assigned to the square for
#'   identification.
#'
#' @return A data frame with columns `Selection`, `X`, and `Y` representing the
#'   square's name and corner coordinates.
#' @export
#'
#' @examples
#' getSquareCoordinates(centerX = 10, centerY = 20, width = 5, name = "spot1")
#' getSquareCoordinates(centerX = 5, centerY = 5, width = 4, name = "MySquare")
getSquareCoordinates <- function(centerX, centerY, width, name) {
    # Calculate half width
    half_width <- width / 2

    # Calculate coordinates of the four corners
    x_coords <- c(
        centerX - half_width,
        centerX + half_width,
        centerX + half_width,
        centerX - half_width,
        centerX - half_width
    )
    y_coords <- c(
        centerY - half_width,
        centerY - half_width,
        centerY + half_width,
        centerY + half_width,
        centerY - half_width
    )

    # Return the coordinates as a matrix
    square_df <- data.frame(
        Selection = rep(name, length(x_coords)),
        X = x_coords,
        Y = y_coords
    )
    return(square_df)
}


rotate <- function(
    angle,
    centerCur
) {
    alpha <- 2 * pi * (angle / 360)
    # center.cur <- c(center.cur, 0)
    # points(center.cur[1], center.cur[2], col = "red")
    tr <- rigid.transl(-centerCur[1], -centerCur[2])
    tr <- rigid.transf(centerCur[1], centerCur[2], alpha) %*% tr
    return(tr)
}


#' Creates a transformation matrix that translates an object
#' in 2D
#'
#' @param translateX,translateY translation of x, y coordinates
#' @noRd

translate <- function(
    translateX,
    translateY
) {
    tr <- rigid.transl(translateX, translateY)
    return(tr)
}


#' Creates a transformation matrix that mirrors an object
#' in 2D along either the x axis or y axis around its
#' center of mass
#'
#' @param mirrorX,mirrorY Logical specifying whether or not an
#' object should be reflected
#' @param centerCur Coordinates of the current center of mass
#'
#' @noRd

mirror <- function(
    mirrorX = FALSE,
    mirrorY = FALSE,
    centerCur
) {
    centerCur <- c(centerCur, 0)
    tr <- rigid.transl(-centerCur[1], -centerCur[2])
    tr <- rigid.refl(mirrorX, mirrorY) %*% tr
    tr <- rigid.transl(centerCur[1], centerCur[2]) %*% tr
    return(tr)
}


#' Stretch along angle
#'
#' Creates a transformation matrix that stretches an object
#' along a specific axis
#'
#' @param r stretching factor
#' @param alpha angle
#' @param centerCur Coordinates of the current center of mass
#'
#' @noRd

stretch <- function(r, alpha, centerCur) {
    centerCur <- c(centerCur, 0)
    tr <- rigid.transl(-centerCur[1], -centerCur[2])
    tr <- rigid.rot(alpha, forward = TRUE) %*% tr
    tr <- rigid.stretch(r) %*% tr
    tr <- rigid.rot(alpha, forward = FALSE) %*% tr
    tr <- rigid.transl(centerCur[1], centerCur[2]) %*% tr
    return(tr)
}


#' Creates a transformation matrix for rotation
#'
#' Creates a transformation matrix for clockwise rotation by 'alpha' degrees
#'
#' @param alpha rotation angle
#' @param forward should the rotation be done in forward direction?
#'
#' @noRd

rigid.rot <- function(
    alpha = 0,
    forward = TRUE
) {
    alpha <- 2 * pi * (alpha / 360)
    tr <- matrix(
        c(
            cos(alpha),
            ifelse(forward, -sin(alpha), sin(alpha)),
            0,
            ifelse(forward, sin(alpha), -sin(alpha)),
            cos(alpha),
            0,
            0,
            0,
            1
        ),
        nrow = 3
    )
    return(tr)
}


#' Creates a transformation matrix for rotation and translation
#'
#' Creates a transformation matrix for clockwise rotation by 'alpha' degrees
#' followed by a translation with an offset of (h, k). Points are assumed to be
#' centered at (0, 0).
#'
#' @param h Numeric: offset along x axis
#' @param k Numeric: offset along y axis
#' @param alpha rotation angle
#'
#' @noRd

rigid.transf <- function(
    h = 0,
    k = 0,
    alpha = 0
) {
    tr <- matrix(
        c(cos(alpha), -sin(alpha), 0, sin(alpha), cos(alpha), 0, h, k, 1),
        nrow = 3
    )
    return(tr)
}

#' Creates a transformation matrix for translation with an offset of (h, k)
#'
#' @param h Numeric: offset along x axis
#' @param k Numeric: offset along y axis
#'
#' @noRd

rigid.transl <- function(
    h = 0,
    k = 0
) {
    tr <- matrix(c(1, 0, 0, 0, 1, 0, h, k, 1), nrow = 3)
    return(tr)
}

#' Creates a transformation matrix for reflection
#'
#' Creates a transformation matrix for reflection where mirror.x will reflect
#' the
#' points along the x axis and mirror.y will reflect thepoints along the y axis.
#' Points are assumed to be centered at (0, 0)
#'
#' @param mirrorX,mirrorY Logical: mirrors x or y axis if set to TRUE
#' @noRd

rigid.refl <- function(
    mirrorX,
    mirrorY
) {
    tr <- diag(c(1, 1, 1))
    if (mirrorX) {
        tr[1, 1] <- -tr[1, 1]
    }
    if (mirrorY) {
        tr[2, 2] <- -tr[2, 2]
    }
    return(tr)
}

#' Creates a transformation matrix for stretching
#'
#' Creates a transformation matrix for stretching by a factor of r
#' along the x axis.
#'
#' @param r stretching factor
#' @noRd

rigid.stretch <- function(
    r
) {
    tr <- matrix(c(r, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
}


#' Combines rigid tranformation matrices
#'
#' Combines rigid tranformation matrices in the following order:
#' translation of points to origin (0, 0) -> reflection of points
#' -> rotation by alpha degrees and translation of points to new center
#'
#' @param centerCur (x, y) image pixel coordinates specifying the current
#'   center of the tissue (stored as the `centers` auxiliary entry)
#' @param centerNew (x, y) image pixel coordinates specifying the new center
#'   (image center)
#' @param alpha Rotation angle
#' @param mirrorX,mirrorY Logical values indicating reflection across the x
#'   or y axis before rotation.
#'
#' @return A 3-by-3 homogeneous affine-transformation matrix.
#'
#' @examples
#' combineTransformations(centerCur = c(0, 0), centerNew = c(10, 20), alpha = 30)
#' transformation <- combineTransformations(
#'     centerCur = c(0, 0),
#'     centerNew = c(10, 20),
#'     alpha = 90
#' )
#' transformation
#'
#' @export

combineTransformations <- function(
    centerCur,
    centerNew,
    alpha,
    mirrorX = FALSE,
    mirrorY = FALSE
) {
    alpha <- 2 * pi * (alpha / 360)
    centerCur <- c(centerCur, 0)
    tr <- rigid.transl(-centerCur[1], -centerCur[2])

    # reflect
    tr <- rigid.refl(mirrorX, mirrorY) %*% tr

    # rotate and translate to new center
    tr <- rigid.transf(centerNew[1], centerNew[2], alpha) %*% tr
    tr
}


generate.map.affine <- function(
    tr,
    forward = FALSE
) {
    if (forward) {
        map.affine <- function(x, y) {
            p <- cbind(x, y)
            xy <- t(solve(tr) %*% t(cbind(p, 1)))
            list(x = xy[, 1], y = xy[, 2])
        }
    } else {
        map.affine <- function(x, y) {
            p <- cbind(x, y)
            xy <- t(tr %*% t(cbind(p, 1)))
            list(x = xy[, 1], y = xy[, 2])
        }
    }
    return(map.affine)
}
