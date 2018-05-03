for i=2:3 % two and three state clusterings
    for j=1:10 % ten iterations to average over
        [i,j]
        markov(1,i,j);
    end
end
