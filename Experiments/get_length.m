function length=get_length(sensor_length,R,H)
    theta0=sensor_length/R/4;
    pos1=[0,R]+[0,H];
    pos2=[(R+H)*sin(theta0),(R+H)*cos(theta0)];
    pos3=[(R+H)*sin(2*theta0),(R+H)*cos(2*theta0)];
    
    length=2*(norm(pos1-pos2)+norm(pos2-pos3));
end