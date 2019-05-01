
s = {}

function docalculation (p) 
    local l = s[3] 
    local xoffset = s[0]
    local yoffset = s[1]
    local zoffset = s[2]
    p:addpoint(l, xoffset, yoffset, zoffset)
end
