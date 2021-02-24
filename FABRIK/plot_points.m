function plot_points(p)
    x=p(:,1);
    y=p(:,2);
    plot(x,y,'-o','MarkerSize',2)
    axis equal
    grid on
end