% point cloud

function X = point_cloud(A)
    N = size(A,1);
    Xin = linspace(-20,20,N);
    Yin = Xin;
    n = 100;
    theta = 0 : pi/100 : 2*pi;
    X = [];
    for i = 1:length(theta)
        R = [cos(theta(i)) -sin(theta(i)) 0 ; sin(theta(i)) cos(theta(i)) 0 ; 0 0 1];
            for j = 1:n
                [y,z] = pinky(Xin,Yin,A);
                C = R*[0,y,z]';
                X(end+1,:) = C;
            end
    end
    scatter3(X(:,1),X(:,2),X(:,3),'r.');
    axis equal;
    xlabel('x (angstrom)');
    ylabel('y (angstrom)');
    zlabel('z (angstrom)');
    camproj('perspective');
    view(58,12);
end


