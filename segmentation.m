clear all;
close all;
Im = imresize(imread('2.jpeg'), 0.1);
r = 1.5; 
ImS   = 6;                
SpS   = 6;                 
segNcut = 1;
segArea = 100;
imshow(Im);
[numRow, numCol,c] = size(Im);   
N = numRow * numCol;
V = reshape(Im, N, c);
W = sparse(N,N);       
WI = reshape(Im, N, 1, c); 
Xs = cat(3, repmat((1:numRow)', 1, numCol), repmat((1:numCol), numRow, 1));
Xs = reshape(Xs, N, 1, 2); 
for icol=1:numCol  
    for irow=1:numRow 
        jcol = (icol - floor(r)) : (icol + floor(r));
        jrow = ((irow - floor(r)) :(irow + floor(r)))';
        jcol = jcol(jcol >= 1 & jcol <= numCol);
        jrow = jrow(jrow >= 1 & jrow <= numRow);
        jN = length(jcol) * length(jrow);
        a = irow + (icol - 1) * numRow;
        b = repmat(jrow, 1, length(jcol)) + repmat((jcol -1) * numRow, length(jrow), 1);
        b = reshape(b, length(jcol) * length(jrow), 1);
        XB = Xs(b, 1, :); 
        XA = repmat(Xs(a, 1, :), length(b), 1);
        DXab = XA - XB;
        DXab = sum(DXab .* DXab, 3);
        constraint = find(sqrt(DXab) <= r);
        b = b(constraint);
        DXab = DXab(constraint);
        IB = WI(b, 1, :);
        IA = repmat(WI(a, 1, :), length(b), 1);
        DifIab = IA - IB;
        DifIab = sum(DifIab .* DifIab, 3);
        W(a, b) = exp(-DifIab / (ImS*ImS)) .* exp(-DXab / (SpS*SpS)); 
    end    
end
Seg = (1:N)';
id = 'ROOT';
N = length(W);
d=sum(W,2);
D = spdiags(d, 0, N, N);
[U,S] = eigs(D-W, D, 2, 'sm'); 
U2 = U(:,2); 
t=mean(U2);
t=fminsearch('Value',t,[],U2,W,D);
segA = find(U2 > t);
segB = find(U2 <= t);
x=(U2 > t);
x=(2*x)-1;
d=diag(D);
k=sum(d(x>0))/sum(d);
x = (U2 > t);
x = (2 * x) - 1;
d = diag(D);
k = sum(d(x > 0)) / sum(d);
b = k / (1 - k);
y = (1 + x) - b * (1 - x);
ncut = (y' * (D - W) * y) / ( y' * D * y );
[SegA  NcutA] = Partition(Seg(segA), W(segA, segA), segNcut, segArea );
[SegB  NcutB] = Partition(Seg(segB), W(segB, segB), segNcut, segArea );
Seg  = [SegA SegB];
Ncut = [NcutA NcutB];
NcutImage  = zeros(size(Im),'uint8');
for k=1:length(Seg)
 [r, c] = ind2sub(size(Im),Seg{k});
 for a=1:length(r)
 NcutImage(r(a),c(a),1:3) = uint8(round(mean(V(Seg{k}, :))));
 end
end
figure;
imshow(NcutImage)