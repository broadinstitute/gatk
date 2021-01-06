cat $1 | grep -v "#" | awk '{print $1 ":" $2 }' | uniq > /tmp/f1.pos
cat $2 | grep -v "#" | awk '{print $1 ":" $2 }' | uniq > /tmp/f2.pos

echo "Positions in File 1:" $(cat /tmp/f1.pos | wc -l)
echo "Positions in File 2:" $(cat /tmp/f2.pos | wc -l)

echo "Positions in common:" $(join /tmp/f1.pos /tmp/f2.pos  | wc -l)

diff /tmp/f1.pos /tmp/f2.pos | grep "<" | cut -c3- > /tmp/only.f1.pos
diff /tmp/f1.pos /tmp/f2.pos | grep ">" | cut -c3- > /tmp/only.f2.pos

echo "Positions ONLY IN File 1:" $(cat /tmp/only.f1.pos | wc -l) " -- See /tmp/only.f1.pos"
echo "Positions ONLY IN File 2:" $(cat /tmp/only.f2.pos | wc -l) " -- See /tmp/only.f2.pos"
