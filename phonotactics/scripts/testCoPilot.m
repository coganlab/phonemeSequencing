# find contours in the thresholded image less than 0.5
# and give them a red border
function process_image(image)
    # find contours in the thresholded image
    contours = find_contours(image, 0.5)
    # draw contours in red
    for contour in contours
        draw_contour!(image, contour, color = RGB(1, 0, 0))
    end
    return image
end