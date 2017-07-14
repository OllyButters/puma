var base = d3.select("#network");

var data, tick;
tick = 0;

var width = 1200;
var height = 800;

var renderer = PIXI.autoDetectRenderer(width, height);
renderer.backgroundColor = 0xffffff;

var chart = document.getElementById('network').appendChild(renderer.view);

var world = new PIXI.Container();

var stage = new PIXI.Container();

var highlight_sprite;

world.addChild(stage);

stage.scale.x = 0.5;
stage.scale.y = 0.5;

PIXI.loader
  .add([
  "dot.png",
  "dot-highlight.png",
  "dot4.png",
  "vline3.png",
  "point-40-1px.png"
  ])
  .load(runSim);

//zoom in/out on mouse scroll
function zoomContainer(e) {
  e.preventDefault();
  var zoom_speed = 0.001;
  stage.scale.x -= e.deltaY * zoom_speed;
  stage.scale.y -= e.deltaY * zoom_speed;
  drawGraph(data);
}

//move (repostion stage) when mouse down and moving
var down = {};

chart.addEventListener("mousedown", function(e) {
  simulation.stop();
  down.x = e.clientX;
  down.y = e.clientY;
  chart.addEventListener("mousemove", positionContainer);
  chart.addEventListener("wheel", zoomContainer);
});

chart.addEventListener("mouseup", function(e) {
  chart.removeEventListener("mousemove", positionContainer);
  chart.removeEventListener("wheel", zoomContainer);
  simulation.restart();
});

function positionContainer(e) {
  console.log(e.clientX);
  console.log(e.clientY);
  var delta_x = e.clientX - down.x;
  var delta_y = e.clientY - down.y;
  stage.position.x += delta_x;
  stage.position.y += delta_y;
  down.x = e.clientX;
  down.y = e.clientY;
  drawGraph(data);
}

//at particular zoom level we want to start displaying names
//only add them if visible, so get current co-ords to check if present in renderer width/height and camera pos
function displayLabel() {
}

function testSetup() {
  sprite = new PIXI.Sprite(PIXI.loader.resources["dot.png"]);
  stage.addChild(sprite);
}

function linkStrength(link) {
  return link.count/2;
}

//position in shells around the centre
var shell_step = 100;
function positionX(node) {
  var centre_diff = Math.round((node.x - (width / 2)) / shell_step); 
}

function positionY(node) {
  var centre_diff = Math.round((node.y - (height / 2)) / shell_step); 
}

var simulation = d3.forceSimulation()
  .force("link", d3.forceLink().id(function(d) { return d.id; }))
  .force("charge", d3.forceManyBody())
  .force("center", d3.forceCenter(width / 2, height / 2));
 
function drawSpriteLine() {
  this.sprite.height = Math.sqrt(Math.pow(this.source.x-this.target.x, 2)+Math.pow(this.source.y-this.target.y, 2));

  //sprite rotation is -1 * angle_from_x_axis
  this.sprite.rotation = -1 * Math.atan2(this.target.x-this.source.x, this.target.y-this.source.y);
  this.sprite.x = this.source.x;
  this.sprite.y = this.source.y;
  return this;
}

function moveSprite() {
  this.sprite.x = this.x;
  this.sprite.y = this.y;
}

function drawGraph(data) {
  //position the links
  data.links.forEach(function(l) {
    l.drawSpriteLine();
  });

  //position the nodes
  data.nodes.forEach(function(n) {
    n.moveSprite();
  });
  renderer.render(world);
  tick++;
}

function displayGraph(error, data) {
  if (error !== null) {
    alert("Error: "+error);
    throw error;
  }

  var nodes, links, tick, max_x, max_y, min_x, min_y, scale_x, scale_y;
  
  nodes = data.nodes;
  links = data.links;
  tick = 0;
  max_x = width;
  max_y = height;
  min_x = 0;
  min_y = 0;
  scale_x = 0;
  scale_y = 0;

  var simDrawGraph = drawGraph.bind(undefined, data);

  simulation
    .nodes(data.nodes)
    .on("tick", simDrawGraph);

  simulation.force("link")
    .links(data.links);

  simulation.force("link")
    .strength(linkStrength);

  simulation.force("charge")
    .strength(-10);

  //pixi sprite setup
  function setupSprites(element, type) {
    var sprite;
    switch (type) {
      case "node":
        sprite = new PIXI.Sprite(PIXI.loader.resources["dot.png"].texture);
        sprite.anchor.x = 0.5;
        sprite.anchor.y = 0.5;
        if (typeof element.affiliation !== "undefined" && element.affiliation.length > 0) { 
          var filter = new PIXI.filters.ColorMatrixFilter();
          sprite.filters = [filter];
          filter.hue(data.affiliations[element.affiliation[0]].colour, true);
        }
      break;
      case "link":
        sprite = new PIXI.Sprite(PIXI.loader.resources["point-40-1px.png"].texture);
        sprite.anchor.x = 0.5;
        sprite.anchor.y = 0;
        sprite.alpha = 0.8;
      break;
    }
    return sprite;
  }

  //setup sprite in pixi container for each node and link
  data.links.forEach(function(l) {
    l.sprite = setupSprites(l, "link");
    l.drawSpriteLine = drawSpriteLine;
    stage.addChild(l.sprite);
  });

  data.nodes.forEach(function(n) {
    n.sprite = setupSprites(n, "node");
    n.moveSprite = moveSprite;
    stage.addChild(n.sprite);
    //get a count of links for each node
    var link_count = d3.sum(data.links.map(function(link) {
      if (link.source.id == n.id || link.target.id == n.id) {
        return link.count;
      }
    }));
    n.link_count = link_count;
  });
  
  addKeyTable(data);
  addTopLinksTable(data);

  highlight_sprite = new PIXI.Sprite(PIXI.loader.resources["dot-highlight.png"].texture);
  stage.addChild(highlight_sprite);
  highlight_sprite.anchor.x = 0.5;
  highlight_sprite.anchor.y = 0.5;
}

//add table for affiliations key
function addKeyTable(data) {
  var key_tab = document.createElement('table');
  var row = key_tab.insertRow();
  var th = document.createElement('th');
  th.appendChild(document.createTextNode('Colour'));
  row.appendChild(th);
  th = document.createElement('th');
  th.appendChild(document.createTextNode('Affiliation'));
  row.appendChild(th);

  for (var affil in data.affiliations) {
    row = key_tab.insertRow();
    var cell = row.insertCell();
    var colour = document.createElement('div');
    cell.style.backgroundColor = '#'+data.affiliations[affil].colour;
    colour.width = '10px';
    colour.height = '10px';
    colour.position = 'absolute';
    cell.appendChild(colour);
    cell = row.insertCell();
    cell.appendChild(document.createTextNode(affil));
  }
  document.getElementById('key').appendChild(key_tab);
}

//add table for top linkers
function addTopLinksTable(data) {
  var toplink_tab = document.createElement('table');
  var row = toplink_tab.insertRow();
  var th = document.createElement('th');
  th.appendChild(document.createTextNode('Name'));
  row.appendChild(th);
  th = document.createElement('th');
  th.appendChild(document.createTextNode('Links'));
  row.appendChild(th);

  data.nodes.sort(function(n1, n2) {
    return n2.link_count - n1.link_count;
  });

  data.nodes.forEach(function(node) {
    row = toplink_tab.insertRow();
    var cell = row.insertCell();
    cell.appendChild(document.createTextNode(node.clean_value));
    cell.dataset.node_id = node.id;
    cell.addEventListener('click', highlightAuthor);
    cell = row.insertCell();
    cell.appendChild(document.createTextNode(node.link_count));
  });
  document.getElementById('top_links').appendChild(toplink_tab);
}

function highlightAuthor(e) {
  var node_id = e.target.dataset.node_id;
  console.log(node_id);
  var node = {};
  data.nodes.forEach(function(n) {
    if (n.id == node_id) {
      node = n;
      return;
    }
  });
  highlight_sprite.x = node.x;
  highlight_sprite.y = node.y;
  renderer.render(world);
}

function runSim(datafile) {
  var data_queue = d3.queue();
  /*
  data_queue.defer(d3.json, 'network-nodes-assoc-subject-meshheading.json')
    .defer(d3.json, 'network-assoc_nodes-assoc-subject-meshheading.json')
    .defer(d3.json, 'network-links-assoc-subject-meshheading.json')
    .await(passData);
  */
  data_queue.defer(d3.json, datafile)
    .await(passData);

  function passData(parse_error, node_data) {
    node_data.nodes = node_data.nodes.map(function(node) {
      node.type = 'node';
      return node;
    });
    node_data.assoc_nodes = node_data.assoc_nodes.map(function(node) {
      node.type = 'assoc_node';
      return node;
    });
    node_data.nodes = node_data.nodes.concat(node_data.assoc_nodes);
    node_data.links = node_data.links.map(function(link) {
      return {
        'source': link.src_id,
        'target': link.target_id,
        'count': link.count,
        'id': link.id
      };
    });
    var affiliations = {};
    for (var n = 0; n < node_data.nodes.length; n++) {
      var affiliation = node_data.nodes[n].affiliation;
      if (affiliation.length > 0) {
        for (var a in affiliation) {
          if (affiliations.hasOwnProperty(affiliation[a]) === false) {
            affiliations[affiliation[a]] = {
              'colour': Math.floor(Math.random()*1000)+50,
              'count': 1,
            };
          } else {
            affiliations[affiliation[a]].count++;
          }
        }
      }
    }
    node_data.affiliations = affiliations;
    var affiliation_links = {};
    console.log(node_data);
    data = node_data
    displayGraph(parse_error, node_data);
  }
}
