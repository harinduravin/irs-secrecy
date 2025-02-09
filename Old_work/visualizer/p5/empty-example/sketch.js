class Circle {
  constructor(x,y,r,text,name) {
    this.p = createVector(x,y)
    this.r = r
    this.text = text
    this.name = name
  }
  
  draw() {
    textSize(15);
    text(this.text, this.p.x-10, this.p.y+5);
    fill(10, 101, 192, 127);
    stroke(127, 63, 120);
    circle(this.p.x, this.p.y, this.r)
  }
  
}

let bg

function setup() {
  hover = null
  grabbed = null
  tx_circle = null
  rx_circle = null
  irs_circle = null
  h = 503
  w = 702
  bg = loadImage('background.png')
  createCanvas(w, h)
  ellipseMode(RADIUS)
  circles = []

  circles.push(new Circle(random(width), random(height), 20,"A1","A1"))
  circles.push(new Circle(random(width), random(height), 20,"A2","A2"))
  circles.push(new Circle(random(width), random(height), 20,"B1","B1"))
  circles.push(new Circle(random(width), random(height), 20,"B2","B2"))
  circles.push(new Circle(random(width), random(height), 20,"IRS","IRS"))
  circles.push(new Circle(random(width), random(height), 20,"C","C"))
}

function draw() {
  background(bg)
  m = createVector(mouseX, mouseY)
  hover = null
  for (let c of circles) {
    if (m.dist(c.p) < c.r) {
      hover = c
    }
  }
  // background('white')
  noStroke()
  if (hover) cursor('grab')
  else cursor(ARROW)
  for (let c of circles) {
    if (c == grabbed) fill(50)
    else if (c == hover) fill(100)
    else fill(0)
    c.draw()
  }

}

function mousePressed() {
  if (hover) {
    grabbed = hover
  }
}

function mouseReleased() {
  console.log('HAA = get_H('+get_value_str('A1','A2')+',lambda,L,Lo,pl,10*seed);\r\nHBB = get_H('+get_value_str('B1','B2')+',lambda,L,Lo,pl,10*seed+1);\r\nH11 = get_H('+get_value_str('A1','B1')+',lambda,L,Lo,pl,10*seed+2);\r\n H21 = get_H('+get_value_str('A2','B1')+',lambda,L,Lo,pl,10*seed+3);\r\nH12 = get_H('+get_value_str('A1','B2')+',lambda,L,Lo,pl,10*seed+4);\r\nH22 = get_H('+get_value_str('A2','B2')+',lambda,L,Lo,pl,10*seed+5);\r\nHA1C = get_H('+get_value_str('A1','C')+',lambda,L,Lo,pl,10*seed+6);\r\nHA2C = get_H('+get_value_str('A2','C')+',lambda,L,Lo,pl,10*seed+7);\r\nHB1C = get_H('+get_value_str('B1','C')+',lambda,L,Lo,pl,10*seed+8);\r\nHB2C = get_H('+get_value_str('B2','C')+',lambda,L,Lo,pl,10*seed+9);')


  // console.log(2*(hover.p.x-350)/25)
  // console.log(2*(250 - hover.p.y)/25)
  grabbed = null
}

function mouseDragged() {
  if (grabbed) {
    grabbed.p.add(createVector(movedX, movedY))
  }
}

function get_value_str(tx,rx) {
  for (let c of circles) {
    if (c.name == tx) {
      tx_circle = c
    }
    if (c.name == rx) {
      rx_circle = c
    }
    if (c.name == "IRS") {
      irs_circle = c
    }
  }
  ret_string = get_distance(tx_circle.p,irs_circle.p)+','+get_distance(irs_circle.p,rx_circle.p)+','+get_distance(tx_circle.p,rx_circle.p)
  return ret_string
}

function get_distance(tx_vec, rx_vec) {
  distance = sqrt(Math.pow(tx_vec.x-rx_vec.x,2)+Math.pow(tx_vec.y-rx_vec.y,2))*2/25

  return distance.toFixed(2)
}
