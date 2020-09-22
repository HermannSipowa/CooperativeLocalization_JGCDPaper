function RN = ReferenceFrame(r_chaser,r_target)
% Computing the reference frame
Nrhovec = r_target - r_chaser ; n1 = [1 0 0]'; n3 = [0 0 1]';
r3 = Nrhovec/norm(Nrhovec);
if norm(cross(r3,n3)) ~= 0
    r1 = cross(r3,n3)/norm(cross(r3,n3));
else
    if dot(r3,n3)<0
        r1 = cross(r3,n1)/norm(cross(r3,-n1));
    else
        r1 = cross(r3,n1)/norm(cross(r3,n1));
    end
end
r2 = cross(r3,r1)/norm(cross(r3,r1));
RN = [r1'; r2';r3'];
end