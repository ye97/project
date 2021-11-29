function p = transform(p, R, t)

p = p*R;

            % translate
p(:,1) = p(:,1) + t(1);
p(:,2) = p(:,2) + t(2);
p(:,3) = p(:,3) + t(3);